/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.filter.CountingDuplicateFilter;
import picard.filter.CountingFilter;
import picard.filter.CountingMapQFilter;
import picard.filter.CountingPairedFilter;
import picard.util.MathUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.logging.Level;
import java.util.logging.Logger;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * @author tfennell
 */
@CommandLineProgramProperties(
        usage = CollectWgsMetrics.USAGE_SUMMARY + CollectWgsMetrics.USAGE_DETAILS,
        usageShort = CollectWgsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectWgsMetrics extends CommandLineProgram {
static final String USAGE_SUMMARY = "Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.";
static final String USAGE_DETAILS = "<p>This tool collects metrics about the fractions of reads that pass base- and mapping-quality "+
"filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum "+
"read depths (coverage cap) are user defined.</p>" +

"<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>"+
"<h4>Usage Example:</h4>"+
"<pre>"  +
"java -jar picard.jar CollectWgsMetrics \\<br /> " +
"      I=input.bam \\<br /> "+
"      O=collect_wgs_metrics.txt \\<br /> " +
"      R=reference_sequence.fasta " +
"</pre>" +
"Please see "+
"<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics'>CollectWgsMetrics</a> "+
"for detailed explanations of the output metrics." +
"<hr />"
;

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence fasta aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage. N bases will be treated as having a base quality " +
            "of negative infinity and will therefore be excluded from coverage regardless of the value of this parameter.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(shortName = "CAP", doc = "Treat positions with coverage exceeding this value as if they had coverage at this value (but calculate the difference for PCT_EXC_CAPPED).", overridable = true)
    public int COVERAGE_CAP = 250;

    @Option(doc="At positions with coverage exceeding this value, completely ignore reads that accumulate beyond this value (so that they will not be considered for PCT_EXC_CAPPED).  Used to keep memory consumption in check, but could create bias if set too low", overridable = true)
    public int LOCUS_ACCUMULATION_CAP = 100000;

    @Option(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Option(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @Option(doc="If true, count unpaired reads, and paired reads with one end unmapped")
    public boolean COUNT_UNPAIRED = false;

    @Option(doc="Sample Size used for Theoretical Het Sensitivity sampling. Default is 10000.", optional = true)
    public int SAMPLE_SIZE=10000;

    @Option(doc = "An interval list file that contains the positions to restrict the assessment. Please note that " +
            "all bases of reads that overlap these intervals will be considered, even if some of those bases extend beyond the boundaries of " +
            "the interval. The ideal use case for this argument is to use it to restrict the calculation to a subset of (whole) contigs. To " +
            "restrict the calculation just to individual positions without overlap, please see CollectWgsMetricsFromSampledSites.",
            optional = true, overridable = true)
    public File INTERVALS = null;

    private SAMFileHeader header = null;

    private final Log log = Log.getInstance(CollectWgsMetrics.class);
    private static final double LOG_ODDS_THRESHOLD = 3.0;

    public static final int READS_IN_PACK = 1000;

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MergeableMetricBase {

        /** The intervals over which this metric was computed. */
        @MergingIsManual
        protected IntervalList intervals;

        /** Count of sites with a given depth of coverage. Excludes bases with quality below MINIMUM_BASE_QUALITY (default 20) */
        @MergingIsManual
        protected final Histogram<Integer> highQualityDepthHistogram;

        /** Count of sites with a given depth of coverage. Includes all but quality 2 bases */
        @MergingIsManual
        protected final Histogram<Integer> unfilteredDepthHistogram;

        /** The count of bases observed with a given base quality. */
        @MergingIsManual
        protected final Histogram<Integer> unfilteredBaseQHistogram;

        /** The maximum depth/coverage to consider. */
        @MergeByAssertEquals
        protected final int coverageCap;

        /** The sample size used for theoretical het sensitivity. */
        @NoMergingKeepsValue
        protected final int theoreticalHetSensitivitySampleSize;

        /**
         * Create an instance of this metric that is not mergeable.
         */
        public WgsMetrics() {
            intervals                           = null;
            highQualityDepthHistogram           = null;
            unfilteredDepthHistogram            = null;
            unfilteredBaseQHistogram            = null;
            theoreticalHetSensitivitySampleSize = -1;
            coverageCap                         = -1;
        }

        /**
         * Create an instance of this metric that is mergeable.
         *
         * @param highQualityDepthHistogram the count of genomic positions observed for each observed depth. Excludes bases with quality below MINIMUM_BASE_QUALITY.
         * @param unfilteredDepthHistogram the depth histogram that includes all but quality 2 bases.
         * @param pctExcludedByMapq the fraction of aligned bases that were filtered out because they were in reads with low mapping quality.
         * @param pctExcludedByDupes the fraction of aligned bases that were filtered out because they were in reads marked as duplicates.
         * @param pctExcludedByPairing the fraction of bases that were filtered out because they were in reads without a mapped mate pair.
         * @param pctExcludedByBaseq the fraction of aligned bases that were filtered out because they were of low base quality.
         * @param pctExcludedByOverlap the fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
         * @param pctExcludedByCapping the fraction of aligned bases that were filtered out because they would have raised coverage above the capped value.
         * @param pctExcludeTotal the fraction of bases excluded across all filters.
         * @param coverageCap Treat positions with coverage exceeding this value as if they had coverage at this value.
         * @param unfilteredBaseQHistogram the count of bases observed with a given quality. Includes all but quality 2 bases.
         * @param theoreticalHetSensitivitySampleSize the sample size used for theoretical het sensitivity sampling.
         */
        public WgsMetrics(final IntervalList intervals,
                          final Histogram<Integer> highQualityDepthHistogram,
                          final Histogram<Integer> unfilteredDepthHistogram,
                          final double pctExcludedByMapq,
                          final double pctExcludedByDupes,
                          final double pctExcludedByPairing,
                          final double pctExcludedByBaseq,
                          final double pctExcludedByOverlap,
                          final double pctExcludedByCapping,
                          final double pctExcludeTotal,
                          final int coverageCap,
                          final Histogram<Integer> unfilteredBaseQHistogram,
                          final int theoreticalHetSensitivitySampleSize) {
            this.intervals      = intervals.uniqued();
            this.highQualityDepthHistogram = highQualityDepthHistogram;
            this.unfilteredDepthHistogram = unfilteredDepthHistogram;
            this.unfilteredBaseQHistogram = unfilteredBaseQHistogram;
            this.coverageCap    = coverageCap;
            this.theoreticalHetSensitivitySampleSize = theoreticalHetSensitivitySampleSize;

            PCT_EXC_MAPQ     = pctExcludedByMapq;
            PCT_EXC_DUPE     = pctExcludedByDupes;
            PCT_EXC_UNPAIRED = pctExcludedByPairing;
            PCT_EXC_BASEQ    = pctExcludedByBaseq;
            PCT_EXC_OVERLAP  = pctExcludedByOverlap;
            PCT_EXC_CAPPED   = pctExcludedByCapping;
            PCT_EXC_TOTAL    = pctExcludeTotal;

            calculateDerivedFields();
        }

        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        @NoMergingIsDerived
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        @NoMergingIsDerived
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        @NoMergingIsDerived
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        @NoMergingIsDerived
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        @NoMergingIsDerived
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        @NoMergingIsDerived
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        @NoMergingIsDerived
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        @NoMergingIsDerived
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        @NoMergingIsDerived
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        @NoMergingIsDerived
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        @NoMergingIsDerived
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        @NoMergingIsDerived
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 1X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_1X;
        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_100X;

        /** The theoretical HET SNP sensitivity. */
        @NoMergingIsDerived
        public double HET_SNP_SENSITIVITY;

        /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
        @NoMergingIsDerived
        public double HET_SNP_Q;

        /**
         * Merges the various PCT_EXC_* metrics.
         * @param other metric to merge into this one.
         *
         * @return result of merging, also known as "this"
         */
        @Override
        public MergeableMetricBase merge(final MergeableMetricBase other) {
            final WgsMetrics otherMetric = (WgsMetrics) other;

            if (highQualityDepthHistogram == null || otherMetric.highQualityDepthHistogram == null ||
                    unfilteredDepthHistogram == null || otherMetric.unfilteredDepthHistogram == null) {
                throw new PicardException("Depth histogram is required when deriving metrics.");
            }

            // Union the intervals over which bases are called.  They should have no overlaps!
            // NB: interval lists are already uniqued.
            final long genomeTerritory = this.intervals.getBaseCount() + otherMetric.intervals.getBaseCount();
            this.intervals.addall(otherMetric.intervals.getIntervals());
            this.intervals = this.intervals.uniqued();
            if (this.intervals.getBaseCount() != genomeTerritory) {
                throw new PicardException("Trying to merge WgsMetrics calculated on intervals that overlap.");
            }

            // NB:
            // Since: PCT_EXC_TOTAL     = (totalWithExcludes - thisMetricTotal) / totalWithExcludes;
            // Thus:  totalWithExcludes = total / (1 - PCT_EXC_TOTAL)
            // Proof: Exercise is left to the reader.
            final long thisMetricTotal        = (long) highQualityDepthHistogram.getSum();
            final long otherMetricTotal       = (long) otherMetric.highQualityDepthHistogram.getSum();
            final long total                  = thisMetricTotal + otherMetricTotal;
            final long thisTotalWithExcludes  = (long) (thisMetricTotal / (1.0 - PCT_EXC_TOTAL));
            final long otherTotalWithExcludes = (long) (otherMetricTotal / (1.0 - otherMetric.PCT_EXC_TOTAL));
            final double totalWithExcludes    = thisTotalWithExcludes + otherTotalWithExcludes;

            if (0 < totalWithExcludes) {
                PCT_EXC_DUPE     = (PCT_EXC_DUPE * thisTotalWithExcludes + otherMetric.PCT_EXC_DUPE * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_MAPQ     = (PCT_EXC_MAPQ * thisTotalWithExcludes + otherMetric.PCT_EXC_MAPQ * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_UNPAIRED = (PCT_EXC_UNPAIRED * thisTotalWithExcludes + otherMetric.PCT_EXC_UNPAIRED * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_BASEQ    = (PCT_EXC_BASEQ * thisTotalWithExcludes + otherMetric.PCT_EXC_BASEQ * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_OVERLAP  = (PCT_EXC_OVERLAP * thisTotalWithExcludes + otherMetric.PCT_EXC_OVERLAP * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_CAPPED   = (PCT_EXC_CAPPED * thisTotalWithExcludes + otherMetric.PCT_EXC_CAPPED * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;
            }

            // do any merging that are dictated by the annotations.
            super.merge(other);

            // merge the histograms
            highQualityDepthHistogram.addHistogram(otherMetric.highQualityDepthHistogram);
            unfilteredDepthHistogram.addHistogram(otherMetric.unfilteredDepthHistogram);
            if (unfilteredBaseQHistogram != null && otherMetric.unfilteredBaseQHistogram != null)
                unfilteredBaseQHistogram.addHistogram(otherMetric.unfilteredBaseQHistogram);
            return this;
        }

        @Override
        public void calculateDerivedFields() {
            if (highQualityDepthHistogram == null || unfilteredDepthHistogram == null) throw new PicardException("Depth histogram is required when deriving metrics.");
            if (unfilteredBaseQHistogram != null && theoreticalHetSensitivitySampleSize <= 0) {
                throw new PicardException("Sample size is required when a baseQ histogram is given when deriving metrics.");
            }

            final long[] depthHistogramArray = new long[coverageCap+1];

            for (final Histogram.Bin<Integer> bin : highQualityDepthHistogram.values()) {
                final int depth = Math.min((int) bin.getIdValue(), coverageCap);
                depthHistogramArray[depth] += bin.getValue();
            }

            GENOME_TERRITORY = (long) highQualityDepthHistogram.getSumOfValues();
            MEAN_COVERAGE    = highQualityDepthHistogram.getMean();
            SD_COVERAGE      = highQualityDepthHistogram.getStandardDeviation();
            MEDIAN_COVERAGE  = highQualityDepthHistogram.getMedian();
            MAD_COVERAGE     = highQualityDepthHistogram.getMedianAbsoluteDeviation();

            PCT_1X    = MathUtil.sum(depthHistogramArray, 1, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
            PCT_5X    = MathUtil.sum(depthHistogramArray, 5, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
            PCT_10X   = MathUtil.sum(depthHistogramArray, 10, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_15X   = MathUtil.sum(depthHistogramArray, 15, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_20X   = MathUtil.sum(depthHistogramArray, 20, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_25X   = MathUtil.sum(depthHistogramArray, 25, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_30X   = MathUtil.sum(depthHistogramArray, 30, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_40X   = MathUtil.sum(depthHistogramArray, 40, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_50X   = MathUtil.sum(depthHistogramArray, 50, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_60X   = MathUtil.sum(depthHistogramArray, 60, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_70X   = MathUtil.sum(depthHistogramArray, 70, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_80X   = MathUtil.sum(depthHistogramArray, 80, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_90X   = MathUtil.sum(depthHistogramArray, 90, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_100X  = MathUtil.sum(depthHistogramArray, 100, depthHistogramArray.length) / (double) GENOME_TERRITORY;

            // Get Theoretical Het SNP Sensitivity
            if (unfilteredBaseQHistogram != null && unfilteredDepthHistogram != null) {
                final double[] depthDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredDepthHistogram);
                final double[] baseQDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredBaseQHistogram);
                HET_SNP_SENSITIVITY = TheoreticalSensitivity.hetSNPSensitivity(depthDoubleArray, baseQDoubleArray, theoreticalHetSensitivitySampleSize, LOG_ODDS_THRESHOLD);
                HET_SNP_Q = QualityUtil.getPhredScoreFromErrorProbability((1 - HET_SNP_SENSITIVITY));
            }
        }
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
    }

    /** Gets the SamReader from which records will be examined.  This will also set the header so that it is available in
     *  */
    protected SamReader getSamReader() {
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        this.header        = in.getFileHeader();
        return in;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
        }

        // it doesn't make sense for the locus accumulation cap to be lower than the coverage cap
        if (LOCUS_ACCUMULATION_CAP < COVERAGE_CAP) {
            log.warn("Setting the LOCUS_ACCUMULATION_CAP to be equal to the COVERAGE_CAP (" + COVERAGE_CAP + ") because it should not be lower");
            LOCUS_ACCUMULATION_CAP = COVERAGE_CAP;
        }

        // Setup all the inputs
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "Processed", "loci");
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = getSamReader();
        final SamLocusIterator iterator = getLocusIterator(in);

        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        // The order in which filters are added matters!
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        if (!COUNT_UNPAIRED) {
            filters.add(pairFilter);
        }
        iterator.setSamFilters(filters);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);
        iterator.setMaxReadsToAccumulatePerLocus(LOCUS_ACCUMULATION_CAP);

        final WgsMetricsCollector collector = getCollector(COVERAGE_CAP, getIntervalsToExamine());

        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        List<SamLocusIterator.LocusInfo> records = new ArrayList<>(READS_IN_PACK);
        ExecutorService service = Executors.newCachedThreadPool();

        // Loop through all the loci
        while (iterator.hasNext()) {
            final SamLocusIterator.LocusInfo info = iterator.next();
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());

            // Check that the reference is not N
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (SequenceUtil.isNoCall(base)) continue;

            records.add(info);
            ++counter;

            // Record progress
            progress.record(info.getSequenceName(), info.getPosition());

            if (records.size() < READS_IN_PACK && iterator.hasNext())
            {
                continue;
            }

            final List<SamLocusIterator.LocusInfo> tmpRecords = records;
            records = new ArrayList<>(READS_IN_PACK);

            // Add read pack to the collector in a separate stream
            service.submit(new Runnable() {
                @Override
                public void run() {
                    for (SamLocusIterator.LocusInfo rec : tmpRecords) {
                        collector.addInfo(rec);
                    }
                }
            });

            // Perhaps stop
            if (usingStopAfter && counter > stopAfter) break;
        }

        service.shutdown();

        try {
            service.awaitTermination(1, TimeUnit.DAYS);

        } catch (InterruptedException ex) {
            Logger.getLogger(SinglePassSamProgram.class.getName()).log(Level.SEVERE, null, ex);
        }

        long unfilteredBaseQHistogramSum = 0;
        long unfilteredDepthHistogramSum = 0;
        for (int i = 0; i < collector.unfilteredBaseQHistogramArray.length(); ++i) {
            unfilteredBaseQHistogramSum += collector.unfilteredBaseQHistogramArray.get(i);
        }
        for (int i = 0; i <= collector.coverageCap; ++i) {
            unfilteredDepthHistogramSum += i*collector.unfilteredDepthHistogramArray.get(i);
        }
        if (unfilteredBaseQHistogramSum != unfilteredDepthHistogramSum) {
            throw new PicardException("updated coverage and baseQ distributions unequally");
        }

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        collector.addToMetricsFile(out, INCLUDE_BQ_HISTOGRAM, dupeFilter, mapqFilter, pairFilter);
        out.write(OUTPUT);

        return 0;
    }

    /** Gets the intervals over which we will calculate metrics. */
    protected IntervalList getIntervalsToExamine() {
        final IntervalList intervals;
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervals = IntervalList.fromFile(INTERVALS);
        } else {
            intervals = new IntervalList(this.header);
            for (final SAMSequenceRecord rec : this.header.getSequenceDictionary().getSequences()) {
                final Interval interval = new Interval(rec.getSequenceName(), 1, rec.getSequenceLength());
                intervals.add(interval);
            }
        }
        return intervals;
    }

    /** This method should only be called after {@link this.getSamReader()} is called. */
    protected SAMFileHeader getSamFileHeader() {
        if (this.header == null) throw new IllegalStateException("getSamFileHeader() was called but this.header is null");
        return this.header;
    }


    protected WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                            final Histogram<Integer> highQualityDepthHistogram,
                                            final Histogram<Integer> unfilteredDepthHistogram,
                                            final double pctExcludedByMapq,
                                            final double pctExcludedByDupes,
                                            final double pctExcludedByPairing,
                                            final double pctExcludedByBaseq,
                                            final double pctExcludedByOverlap,
                                            final double pctExcludedByCapping,
                                            final double pctTotal,
                                            final int coverageCap,
                                            final Histogram<Integer> unfilteredBaseQHistogram,
                                            final int theoreticalHetSensitivitySampleSize) {
        return new WgsMetrics(
                intervals,
                highQualityDepthHistogram,
                unfilteredDepthHistogram,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                unfilteredBaseQHistogram,
                theoreticalHetSensitivitySampleSize
        );
    }

    private WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                          final Histogram<Integer> highQualityDepthHistogram,
                                          final Histogram<Integer> unfilteredDepthHistogram,
                                          final long basesExcludedByMapq,
                                          final long basesExcludedByDupes,
                                          final long basesExcludedByPairing,
                                          final AtomicLong basesExcludedByBaseq,
                                          final AtomicLong basesExcludedByOverlap,
                                          final AtomicLong basesExcludedByCapping,
                                          final int coverageCap,
                                          final Histogram<Integer> unfilteredBaseQHistogram,
                                          final int theoreticalHetSensitivitySampleSize) {
        final double total = highQualityDepthHistogram.getSum();
        final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq.get() + basesExcludedByOverlap.get() + basesExcludedByCapping.get();

        final double pctExcludedByMapq    = totalWithExcludes == 0 ? 0 : basesExcludedByMapq         / totalWithExcludes;
        final double pctExcludedByDupes   = totalWithExcludes == 0 ? 0 : basesExcludedByDupes        / totalWithExcludes;
        final double pctExcludedByPairing = totalWithExcludes == 0 ? 0 : basesExcludedByPairing      / totalWithExcludes;
        final double pctExcludedByBaseq   = totalWithExcludes == 0 ? 0 : basesExcludedByBaseq.get()        / totalWithExcludes;
        final double pctExcludedByOverlap = totalWithExcludes == 0 ? 0 : basesExcludedByOverlap.get()      / totalWithExcludes;
        final double pctExcludedByCapping = totalWithExcludes == 0 ? 0 : basesExcludedByCapping.get()      / totalWithExcludes;
        final double pctTotal             = totalWithExcludes == 0 ? 0 : (totalWithExcludes - total) / totalWithExcludes;

        return generateWgsMetrics(
                intervals,
                highQualityDepthHistogram,
                unfilteredDepthHistogram,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                unfilteredBaseQHistogram,
                theoreticalHetSensitivitySampleSize
        );
    }


    /**
     * If INTERVALS is specified, this will count bases beyond the interval list when the read overlaps the intervals and extends beyond the
     * edge. Ideally INTERVALS should only include regions that have hard edges without reads that could extend beyond the boundary (such as a whole contig).
     */
    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredBases();
    }

    protected SamLocusIterator getLocusIterator(final SamReader in) {
        return (INTERVALS != null) ? new SamLocusIterator(in, IntervalList.fromFile(INTERVALS)) : new SamLocusIterator(in);
    }

    /**
     * @param coverageCap the maximum depth/coverage to consider.
     * @param intervals the intervals over which metrics are collected.
     * @return
     */
    protected WgsMetricsCollector getCollector(final int coverageCap, final IntervalList intervals) {
        return new WgsMetricsCollector(coverageCap, intervals);
    }

    protected class WgsMetricsCollector {

        // Count of sites with a given depth of coverage. Includes all but quality 2 bases.
        // We draw depths from this histogram when we calculate the theoretical het sensitivity.
        protected final AtomicLongArray unfilteredDepthHistogramArray;

        // Count of bases observed with a given base quality. Includes all but quality 2 bases.
        // We draw bases from this histogram when we calculate the theoretical het sensitivity.
        protected final AtomicLongArray unfilteredBaseQHistogramArray;

        // Count of sites with a given depth of coverage. Excludes bases with quality below MINIMUM_BASE_QUALITY (default 20).
        protected final AtomicLongArray highQualityDepthHistogramArray;

        private AtomicLong basesExcludedByBaseq = new AtomicLong();
        private AtomicLong basesExcludedByOverlap = new AtomicLong();
        private AtomicLong basesExcludedByCapping = new AtomicLong();
        protected final IntervalList intervals;
        protected final int coverageCap;

        public WgsMetricsCollector(final int coverageCap, final IntervalList intervals) {
            unfilteredDepthHistogramArray = new AtomicLongArray(coverageCap + 1);
            highQualityDepthHistogramArray = new AtomicLongArray(coverageCap + 1);
            unfilteredBaseQHistogramArray = new AtomicLongArray(coverageCap + 1);
            this.coverageCap    = coverageCap;
            this.intervals      = intervals;
        }

        public void addInfo(final SamLocusIterator.LocusInfo info) {

            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            int unfilteredDepth = 0;

            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {
                if (recs.getBaseQuality() <= 2) { basesExcludedByBaseq.incrementAndGet();;   continue; }

                // we add to the base quality histogram any bases that have quality > 2
                // the raw depth may exceed the coverageCap before the high-quality depth does. So stop counting once we reach the coverage cap.
                if (unfilteredDepth < coverageCap) {
                    unfilteredBaseQHistogramArray.getAndIncrement(recs.getRecord().getBaseQualities()[recs.getOffset()]);
                    unfilteredDepth++;
                }

                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY ||
                        SequenceUtil.isNoCall(recs.getReadBase()))                  { basesExcludedByBaseq.incrementAndGet();;   continue; }
                if (!readNames.add(recs.getRecord().getReadName()))                 { basesExcludedByOverlap.incrementAndGet();; continue; }

                pileupSize++;
            }

            final int highQualityDepth = Math.min(pileupSize, coverageCap);
            if (highQualityDepth < pileupSize) basesExcludedByCapping.addAndGet(pileupSize - coverageCap);
            highQualityDepthHistogramArray.getAndIncrement(highQualityDepth);
            unfilteredDepthHistogramArray.getAndIncrement(unfilteredDepth);
        }

        public void addToMetricsFile(final MetricsFile<WgsMetrics, Integer> file,
                                     final boolean includeBQHistogram,
                                     final CountingFilter dupeFilter,
                                     final CountingFilter mapqFilter,
                                     final CountingPairedFilter pairFilter) {
            final WgsMetrics metrics = getMetrics(dupeFilter, mapqFilter, pairFilter);

            // add them to the file
            file.addMetric(metrics);
            file.addHistogram(getHighQualityDepthHistogram());
            if (includeBQHistogram) addBaseQHistogram(file);


        }

        protected void addBaseQHistogram(final MetricsFile<WgsMetrics, Integer> file) {
            file.addHistogram(getUnfilteredBaseQHistogram());
        }

        protected Histogram<Integer> getHighQualityDepthHistogram() {
            return getHistogram(highQualityDepthHistogramArray, "coverage", "high_quality_coverage_count");
        }

        protected Histogram<Integer> getUnfilteredDepthHistogram(){
            return getHistogram(unfilteredDepthHistogramArray, "coverage", "unfiltered_coverage_count");
        }

        protected Histogram<Integer> getUnfilteredBaseQHistogram() {
            return getHistogram(unfilteredBaseQHistogramArray, "baseq", "unfiltered_baseq_count");
        }

        protected Histogram<Integer> getHistogram(final AtomicLongArray array, final String binLabel, final String valueLabel) {
            final Histogram<Integer> histogram = new Histogram<>(binLabel, valueLabel);
            for (int i = 0; i < array.length(); ++i) {
                histogram.increment(i, array.get(i));
            }
            return histogram;
        }

        protected WgsMetrics getMetrics(final CountingFilter dupeFilter,
                                        final CountingFilter mapqFilter,
                                        final CountingPairedFilter pairFilter) {
            return generateWgsMetrics(
                    this.intervals,
                    getHighQualityDepthHistogram(),
                    getUnfilteredDepthHistogram(),
                    getBasesExcludedBy(mapqFilter),
                    getBasesExcludedBy(dupeFilter),
                    getBasesExcludedBy(pairFilter),
                    basesExcludedByBaseq,
                    basesExcludedByOverlap,
                    basesExcludedByCapping,
                    coverageCap,
                    getUnfilteredBaseQHistogram(),
                    SAMPLE_SIZE
            );
        }


    }
}

