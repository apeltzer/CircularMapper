import com.google.common.io.Files;
import htsjdk.samtools.*;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import static org.junit.Assert.assertEquals;


/**
 * Created by peltzer on 13/11/2016.
 */
public class ReadTrimTest {
    private File inputFile;
    private File referencePath;


    /**
     * Setup the test environment
     * @param filepath
     * @param referencepath
     */

    public void setUpTests(String filepath, String referencepath) {
        inputFile = new File(filepath);
        referencePath = new File(referencepath);
    }

    /**
     * Set up the test for a fully in the elongation read
     * @throws IOException
     * @throws CloneNotSupportedException
     */
    @Test
    public void realigner_in_elongation() throws IOException, CloneNotSupportedException {
        setUpTests("src/test/resources/Proper_Subsample_Template_in_elongation.sam", "src/test/resources/hg19_MT.fasta");
        RealignSAMFile.main(new String[]{"-i", inputFile.getAbsolutePath(), "-r", referencePath.getAbsolutePath(), "-e", "500", "-f", "true"});
        String realPath = Files.getNameWithoutExtension(inputFile.getAbsolutePath()) + "_realigned.bam";
        String directory = "src/test/resources/";
        SamReader changedReader = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(directory + "/" + realPath));

        ArrayList<SAMRecord> listReads = new ArrayList<SAMRecord>();

        Iterator it = changedReader.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            listReads.add(curr);
        }


        //Check our Read Number #1

        SAMRecord input_1 = listReads.get(0);

        assertEquals(input_1.getAlignmentStart(), 70);
        assertEquals(input_1.getCigarString(), "65M");
        assertEquals(input_1.getAlignmentEnd(),134);


        //Check our Read Number #2

        SAMRecord input_2 = listReads.get(1);

        assertEquals(input_2.getAlignmentStart(), 1);
        assertEquals(input_2.getCigarString(), "66M");
        assertEquals(input_2.getAlignmentEnd(),66);
    }

    /**
     * Set up the test for a fully in the start read
     */

    @Test
    public void realigner_in_start() throws IOException, CloneNotSupportedException {
        setUpTests("src/test/resources/Proper_Subsample_Template_in_first_n.sam", "src/test/resources/hg19_MT.fasta");
        RealignSAMFile.main(new String[]{"-i", inputFile.getAbsolutePath(), "-r", referencePath.getAbsolutePath(), "-e", "500", "-f", "true"});
        String realPath = Files.getNameWithoutExtension(inputFile.getAbsolutePath()) + "_realigned.bam";
        String directory = "src/test/resources/";
        SamReader changedReader = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(directory + "/" + realPath));

        ArrayList<SAMRecord> listReads = new ArrayList<SAMRecord>();

        Iterator it = changedReader.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            listReads.add(curr);
        }

        //Check our Read Number #1

        SAMRecord input_1 = listReads.get(0);

        assertEquals(input_1.getAlignmentStart(), 1);
        assertEquals(input_1.getCigarString(), "66M");
        assertEquals(input_1.getAlignmentEnd(),66);


        //Check our Read Number #2

        SAMRecord input_2 = listReads.get(1);

        assertEquals(input_2.getAlignmentStart(), 72);
        assertEquals(input_2.getCigarString(), "136M");
        assertEquals(input_2.getAlignmentEnd(),207);

    }

    /**
     * Set up the test for a read that needs to be splitted/soft-clipped correctly.
     */

    @Test
    public void realigner_need_split() throws IOException, CloneNotSupportedException {
        setUpTests("src/test/resources/Proper_Subsample_Template_overlap.sam", "src/test/resources/hg19_MT.fasta");
        RealignSAMFile.main(new String[]{"-i", inputFile.getAbsolutePath(), "-r", referencePath.getAbsolutePath(), "-e", "500", "-f", "true"});
        String realPath = Files.getNameWithoutExtension(inputFile.getAbsolutePath()) + "_realigned.bam";
        String directory = "src/test/resources/";
        SamReader changedReader = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(directory + "/" + realPath));

        ArrayList<SAMRecord> listReads = new ArrayList<SAMRecord>();

        Iterator it = changedReader.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            listReads.add(curr);
        }

        //Check our Read Number #1

        SAMRecord input_1 = listReads.get(0);

        assertEquals(input_1.getAlignmentStart(), 1);
        assertEquals(input_1.getCigarString(), "27S36M");
        assertEquals(input_1.getAlignmentEnd(),36);


        //Check our Read Number #2 (brother of #1)

        SAMRecord input_2 = listReads.get(1);

        assertEquals(input_2.getAlignmentStart(), 16543);
        assertEquals(input_2.getCigarString(), "27M36S");
        assertEquals(input_2.getAlignmentEnd(),16569);

        SAMRecord input_3 = listReads.get(2);

        assertEquals(input_3.getAlignmentStart(), 1);
        assertEquals(input_3.getCigarString(), "20S56M");
        assertEquals(input_3.getAlignmentEnd(),56);

        SAMRecord input_4 = listReads.get(3);

        assertEquals(input_4.getAlignmentStart(), 16550);
        assertEquals(input_4.getCigarString(), "20M56S");
        assertEquals(input_4.getAlignmentEnd(),16569);
    }


    @Test
    public void realigner_in_regular_genome() throws IOException, CloneNotSupportedException {
        setUpTests("src/test/resources/Proper_Subsample_Template_in_regular_genome.sam", "src/test/resources/hg19_MT.fasta");
        RealignSAMFile.main(new String[]{"-i", inputFile.getAbsolutePath(), "-r", referencePath.getAbsolutePath(), "-e", "500", "-f", "true"});
        String realPath = Files.getNameWithoutExtension(inputFile.getAbsolutePath()) + "_realigned.bam";
        String directory = "src/test/resources/";
        SamReader changedReader = SamReaderFactory.make().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(directory + "/" + realPath));

        ArrayList<SAMRecord> listReads = new ArrayList<SAMRecord>();

        Iterator it = changedReader.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            listReads.add(curr);
        }

        //Check our Read Number #1

        SAMRecord input_1 = listReads.get(0);

        assertEquals(input_1.getAlignmentStart(), 882);
        assertEquals(input_1.getCigarString(), "55M");
        assertEquals(input_1.getAlignmentEnd(),936);


        //Check our Read Number #2

        SAMRecord input_2 = listReads.get(1);

        assertEquals(input_2.getAlignmentStart(), 859);
        assertEquals(input_2.getCigarString(), "58M");
        assertEquals(input_2.getAlignmentEnd(),916);

    }




}

