package picard.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

public class VcfTestUtils {

    /**
     * This method creates a temporary VCF file and it's appropriately named index file, and will delete them on exit.
     * @param prefix - The prefix string to be used in generating the file's name; must be at least three characters long
     * @param suffix - The suffix string to be used in generating the file's name; may be null, in which case the suffix ".tmp" will be used
     * @return A File object referencing the newly created temporary VCF file
     * @throws IOException - if a file could not be created.
     */
    public static File createTemporaryIndexedVcfFile(final String prefix, final String suffix) throws IOException {
        final File out = File.createTempFile(prefix, suffix);
        out.deleteOnExit();
        String indexFileExtension = null;
        if (suffix.endsWith("vcf.gz")) {
            indexFileExtension = ".tbi";
        }
        else if (suffix.endsWith("vcf")) {
            indexFileExtension = ".idx";
        }
        if (indexFileExtension != null) {
            final File indexOut = new File(out.getAbsolutePath() + indexFileExtension);
            indexOut.deleteOnExit();
        }
        return out;
    }

    /**
     * Useful test method.  Creates a (temporary) indexed VCF so that we don't have to store the index file in the testdata set.
     * @param vcfFile the vcf file to index
     * @return File a vcf file (index file is created in same path).
     */
    public static File createIndexedVcf(final File vcfFile, final String tempFilePrefix) throws IOException {
        final File output = File.createTempFile(tempFilePrefix, ".vcf");
        output.deleteOnExit();
        final File indexFile = new File(output.getAbsolutePath() + ".idx");
        indexFile.deleteOnExit();
        final VCFFileReader in = new VCFFileReader(vcfFile, false);
        final VCFHeader header = in.getFileHeader();

        final VariantContextWriter out = new VariantContextWriterBuilder().
                setReferenceDictionary(header.getSequenceDictionary()).
                setOptions(EnumSet.of(Options.INDEX_ON_THE_FLY)).
                setOutputFile(output).build();
        out.writeHeader(header);
        for (final VariantContext ctx : in) {
            out.add(ctx);
        }
        out.close();
        in.close();
        return output;
    }
}
