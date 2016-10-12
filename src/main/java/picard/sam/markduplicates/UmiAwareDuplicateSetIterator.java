/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * This acts as an iterator over duplicate sets.  If a particular duplicate
 * set consists of records that contain UMIs this iterator breaks up a single
 * duplicate set into multiple duplicate based on the content of the UMIs.
 * Since there may be sequencing errors in the UMIs, this class allows for
 * simple error correction based on edit distances between the UMIs.
 *
 * @author fleharty
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Histogram;
import picard.PicardException;
import picard.sam.UmiMetrics;

import java.util.*;

/**
 * UmiAwareDuplicateSetIterator is an iterator that wraps a duplicate set iterator
 * in such a way that each duplicate set may be broken up into subsets according
 * to UMIs in the records.  Some tolerance for errors in the UMIs is allowed, and
 * the degree of this is controlled by the maxEditDistanceToJoin parameter.
 */
class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private final DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private final int maxEditDistanceToJoin;
    private final String umiTag;
    private final String inferredUmiTag;
    private final boolean allowMissingUmis;
    private boolean isOpen = false;
    private UmiMetrics metrics;
    private boolean haveWeSeenFirstRead = false;
    private long duplicateSetsWithUmi = 0;
    private long duplicateSetsWithoutUmi = 0;

    private List<Histogram<String>> observedUmisOfVariedLength = new ArrayList<>();
    private Histogram<String> observedUmis = new Histogram<>();
    private Histogram<String> inferredUmis = new Histogram<>();
    private Histogram<String> umiObservationByDuplicateSet = new Histogram<>();
    private Histogram<Double> brokenByUmiDistribution = new Histogram<>("DuplicateSetsFromUmi", "numberOfDuplicateSets");

    /**
     * Creates a UMI aware duplicate set iterator
     *
     * @param wrappedIterator UMI aware duplicate set iterator is a wrapper
     * @param maxEditDistanceToJoin The edit distance between UMIs that will be used to union UMIs into groups
     * @param umiTag The tag used in the bam file that designates the UMI
     * @param assignedUmiTag The tag in the bam file that designates the assigned UMI
     */

    // Metrics is passed in by reference so that iterator can update it.
    UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int maxEditDistanceToJoin,
                                 final String umiTag, final String assignedUmiTag, final boolean allowMissingUmis,
                                 final UmiMetrics metrics) {
        this.wrappedIterator = wrappedIterator;
        this.maxEditDistanceToJoin = maxEditDistanceToJoin;
        this.umiTag = umiTag;
        this.inferredUmiTag = assignedUmiTag;
        this.allowMissingUmis = allowMissingUmis;
        this.metrics = metrics;
        isOpen = true;
        nextSetsIterator = Collections.emptyIterator();

    }

    @Override
    public void close() {
        isOpen = false;
        wrappedIterator.close();

        metrics.OBSERVED_UNIQUE_UMIS = observedUmis.size();
        metrics.INFERRED_UNIQUE_UMIS = inferredUmis.size();

        metrics.EFFECTIVE_LENGTH_OF_OBSERVED_UMIS = effectiveNumberOfBases(observedUmis);
        metrics.EFFECTIVE_LENGTH_OF_INFERRED_UMIS = effectiveNumberOfBases(inferredUmis);

        metrics.DUPLICATE_SETS_WITH_UMI = duplicateSetsWithUmi;
        metrics.DUPLICATE_SETS_WITHOUT_UMI = duplicateSetsWithoutUmi;

        Histogram<Double> effectiveUmiLengthDistribution = new Histogram<>();

        for(int k = 0; k < metrics.UMI_LENGTH; k++) {
            effectiveUmiLengthDistribution.increment(effectiveNumberOfBases(observedUmisOfVariedLength.get(k)));
        }
        metrics.computeMetrics();
        metrics.setDuplicateSetsBrokenByUmi(brokenByUmiDistribution);
        metrics.setEffectiveUmiLengthDistribution(effectiveUmiLengthDistribution);

    }

    @Override
    public boolean hasNext() {
        if (!isOpen) {
            return false;
        } else {
            if (nextSetsIterator.hasNext() || wrappedIterator.hasNext()) {
                return true;
            } else {
                isOpen = false;
                return false;
            }
        }
    }

    @Override
    public DuplicateSet next() {
        if (!nextSetsIterator.hasNext()) {
            process(wrappedIterator.next());
        }
        return nextSetsIterator.next();
    }

    /**
     * Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
     * and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
     *
     * @param set Duplicate set that may be broken up into subsets according the UMIs
     */
    private void process(final DuplicateSet set) {

        // Ensure that the nextSetsIterator isn't already occupied
        if (nextSetsIterator.hasNext()) {
            throw new PicardException("nextSetsIterator is expected to be empty, but already contains data.");
        }

        final UmiGraph umiGraph = new UmiGraph(set, umiTag, inferredUmiTag, allowMissingUmis);

        List<DuplicateSet> duplicateSets = umiGraph.joinUmisIntoDuplicateSets(maxEditDistanceToJoin);

        for (DuplicateSet ds : duplicateSets) {
            List<SAMRecord> records = ds.getRecords();
            SAMRecord representativeRead = ds.getRepresentative();

            String inferredUmi = representativeRead.getStringAttribute(inferredUmiTag);

            for (SAMRecord rec : records) {
                String currentUmi = rec.getStringAttribute(umiTag);

                if(currentUmi != null) {
                    if (!haveWeSeenFirstRead) {
                        metrics.UMI_LENGTH = currentUmi.length();

                        for (int k = 1; k <= metrics.UMI_LENGTH; k++) {
                            observedUmisOfVariedLength.add(new Histogram<>());
                        }
                        haveWeSeenFirstRead = true;
                    }

                    metrics.OBSERVED_BASE_ERRORS += hammingDistance(currentUmi, inferredUmi);
                    metrics.TOTAL_UMI_BASES_OBSERVED += metrics.UMI_LENGTH;

                    observedUmis.increment(currentUmi);
                    inferredUmis.increment(inferredUmi);

                    for (int k = 0; k < currentUmi.length(); k++) {
                        observedUmisOfVariedLength.get(k).increment(currentUmi.substring(currentUmi.length() - k - 1));
                    }
                }
            }

            // Increment the number of times we have seen a particular barcode
            // in a duplicate set.
            if(inferredUmi != null) {
                umiObservationByDuplicateSet.increment(inferredUmi);
            }

//            if(ds.size() > 100) {
//                measureBaseConcordance(ds);
//            }
            duplicateSetsWithUmi++;
        }
        duplicateSetsWithoutUmi++;

        brokenByUmiDistribution.increment((double) duplicateSets.size());

        nextSetsIterator = duplicateSets.iterator();
    }

    private double effectiveNumberOfBases(Histogram<?> observations) {
        double H = 0.0;

        double totalObservations = observations.getSumOfValues();
        for (Histogram.Bin observation : observations.values()) {
            double p = observation.getValue() / totalObservations;
            H = H - p*Math.log(p);
        }

        // Convert to log base 4 so that the entropy is now a measure
        // of the effective number of DNA bases.  If we used log(2.0)
        // our result would be in bits.
        return H/Math.log(4.0);
    }

    // Calculate Hamming distance between two strings
    // TODO: use HTSJDK version when this become available
    private int hammingDistance(final String s1, final String s2) {
        // Comparing edit distance of strings with different lengths is not supported
        if (s1.length() != s2.length()) {
            throw new PicardException("Attempting to determine if two UMIs of different length were within a specified edit distance.");
        }
        int measuredDistance = 0;
        for (int i = 0;i < s1.length();i++) {
            if (s1.charAt(i) != s2.charAt(i)) {
                measuredDistance++;
            }
        }
        return measuredDistance;
    }

    private double measureBaseConcordance(DuplicateSet ds) {
        List<SAMRecord> records = ds.getRecords();

        byte[][] baseArray = new byte[records.size()][];
        int i = 0;
        for(SAMRecord rec : records) {
            baseArray[i] = rec.getReadBases();
//            System.out.println(new String(baseArray[i]));
            i++;
        }

        for(int j = 0;j < baseArray[0].length;j++) {
            int a = 0;
            int t = 0;
            int c = 0;
            int g = 0;
            for(int k = 0;k < records.size();k++) {
                if(baseArray[k][j] == 'A') {
                    a++;
                }
                if(baseArray[k][j] == 'T') {
                    t++;
                }
                if(baseArray[k][j] == 'C') {
                    c++;
                }
                if(baseArray[k][j] == 'G') {
                    g++;
                }
            }
            System.out.println("A=" + a + "   T=" + t + "   C=" + c + "   G=" + g);
        }
        return 0.0;
    }
}

