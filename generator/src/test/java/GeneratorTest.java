import com.google.common.annotations.VisibleForTesting;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static junit.framework.TestCase.fail;
import static org.junit.Assert.assertEquals;

/**
 * Created by peltzer on 11/11/2016.
 */
public class GeneratorTest {
    File input = new File(getClass().getResource("/ParatyphiC.fasta").getFile());
    File output = null;
    StringBuffer buffer_output = new StringBuffer();
    StringBuffer buffer_input = new StringBuffer();


    public void testSetup() throws IOException {
        //Set things up
        CircularGenerator cg = new CircularGenerator(50);
        cg.setKeys_to_treat_circular("gi|224581838|ref|NC_01");
        cg.extendFastA(input);


        //Read output
        output = new File(getClass().getResource("/ParatyphiC_50.fasta").getFile());
        FileReader fr = new FileReader(output);
        BufferedReader bfr = new BufferedReader(fr);

        //Read input
        FileReader fri = new FileReader(input);
        BufferedReader bfri = new BufferedReader(fri);


        String currline = "";

        while((currline = bfr.readLine()) != null){
            if(!currline.startsWith(">")) {
                buffer_output.append(currline);
            } else {
                continue;
            }
        }

        String currlinei = "";

        while((currlinei = bfri.readLine()) != null){
            if(!currlinei.startsWith(">")) {
                buffer_input.append(currlinei);
            } else {
                continue;
            }
        }
    }


    public void testSetupWrong() throws IOException {
        //Set things up
        CircularGenerator cg = new CircularGenerator(1000000000); //This is much longer than our input sequence and should thus fail!
        cg.setKeys_to_treat_circular("gi|224581838|ref|NC_01");
        cg.extendFastA(input);
    }

    /**
     * This method tests whether our output start is equal to our input start sequence, theye should not have been changed at all.
     * @throws IOException
     */
    @Test
    public void generator_test_start_equal() throws IOException {
        testSetup();
        String first_50_bases_input = buffer_input.substring(0,50);
        String first_50_bases_output = buffer_output.substring(0,50);
        assertEquals(first_50_bases_input,first_50_bases_output);

    }

    /**
     * This method should test whether the elongation works fine. Thus the first 50 bases should be equal to the last 50 bases of our FastA sequence.
     * @throws IOException
     */
    @Test
    public void generator_extension_correct() throws IOException {
        testSetup();
        String first_50_bases_input = buffer_input.substring(0,50);
        int output_length = buffer_output.length();
        String last_50_bases_output = buffer_output.substring(output_length-50,output_length);
        assertEquals(first_50_bases_input,last_50_bases_output);
    }

    @Test (expected = RuntimeException.class)
    public void generator_too_long_extension(){
        try {
            CircularGenerator cg = new CircularGenerator(1000000000); //This is much longer than our input sequence and should thus fail!
            cg.setKeys_to_treat_circular("gi|224581838|ref|NC_01");
            cg.extendFastA(input);
        } catch (RuntimeException re) {
            String message = "You cannot extend your sequence with a value longer than your actual sequence.";
            assertEquals(message, re.getMessage());
            throw re;
        }
        fail("Too long sequence and still got extended by longer than actual sequence.");



    }


}
