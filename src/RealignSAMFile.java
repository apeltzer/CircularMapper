/*
 * Copyright (c) 2016. CircularMapper Alexander Peltzer
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import org.apache.commons.cli.*;
import htsjdk.samtools.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.*;

/**
 * This class is used to realign the shifted SAM files to a given reference
 * genome. The main method requires three parameters: - input file (STRING,
 * filepath) - the factor by which the elongation has been done (INT) - the
 * length (AA bases) of the unmodified reference genome. (INT)
 *
 * @author peltzer
 */
public class RealignSAMFile {
	private static final String CLASS_NAME = "RealignSAMFile";
    private static final String VERSION = "1.0";

    private LinkedList<SAMRecord> mappedReads;
    private final SamReader inputSam;
    private SAMFileHeader modSamheader;
    private SAMFileWriter output_realigned;
    private int elongationfactor;
    private int stat_infactor;
    private int stat_inelongation;
    private int stat_inregulargenome;
    private int stat_overlapping;
    private static boolean filtered = false;
    private static boolean filter_sequence_ids = false;
    private Set<String> changedSequences;
    private ReadTrimmer rt = new ReadTrimmer();


    public RealignSAMFile(File f, File reference, int elongationfactor, String filter) {
        mappedReads = new LinkedList<SAMRecord>();
        inputSam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(f);
        this.elongationfactor = elongationfactor;

        File outreal = new File(f.getAbsolutePath().substring(0, f.getAbsolutePath().length() - 4) + "_realigned.bam");
        File changedSequences = new File(reference.getAbsoluteFile()+"_"+this.elongationfactor+"_elongated");
        // get the sequences that were elongated
        getChangedSequences(changedSequences);
        
        SAMSequenceDictionary dict = modifyHeader();

        //Set the corrected header in the output file
        this.modSamheader = inputSam.getFileHeader().clone();
        this.modSamheader.setSequenceDictionary(dict);

        //Write output to realigned_file
        this.output_realigned = new SAMFileWriterFactory().makeSAMOrBAMWriter(modSamheader, false, outreal);


    }
    
    private SAMSequenceDictionary modifyHeader(){
    	SAMSequenceDictionary dict = inputSam.getFileHeader().getSequenceDictionary();
        SAMSequenceDictionary newDict = new SAMSequenceDictionary();
    	for(SAMSequenceRecord seq: dict.getSequences()){
    		if(this.changedSequences.contains(seq.getSequenceName())){
    			SAMSequenceRecord rec = seq.clone();
    			int unmodLength = rec.getSequenceLength();
    			int modLength = unmodLength - this.elongationfactor;
    			rec.setSequenceLength(modLength);
    			newDict.addSequence(rec);
                //lis.add(rec);
    		} else if(!filter_sequence_ids){
                newDict.addSequence(seq.clone());
               // lis.add(seq);
            } else {
            }
    	}
    	return newDict;
    }

	private void getChangedSequences(File changedSequences) {
    	this.changedSequences = new HashSet<String>();
    	try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(changedSequences));
			String currLine = "";
			while((currLine = br.readLine()) != null){
				this.changedSequences.add(currLine);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("static-access")
	public static void main(String[] args) throws IOException, CloneNotSupportedException {

		RealignSAMFile rs =  addCLIInterface(args);
		rs.readSAMFile();
		rs.output_realigned.close();
		rs.inputSam.close();
        /**
         * Some statistics output
         */

        System.out.println("Number of overlapping reads: " + rs.stat_overlapping + "\n");
        System.out.println("Number of reads in regular genome: " + rs.stat_inregulargenome + "\n");
        System.out.println("Number of reads in elongation factor completely: " + rs.stat_inelongation + "\n");
        System.out.println("Number of reads in first n bases, where n is the length of the elongation: " + rs.stat_infactor + "\n");

    }



    /**
     * This method reads a SAM file.
     */
    private void readSAMFile() throws IOException, CloneNotSupportedException {
        Iterator it = inputSam.iterator();
        while (it.hasNext()) {
            SAMRecord curr = (SAMRecord) it.next();
            inspectReads(curr);
        }

    }

    /**
     * This method trims a certain read according to its relative position in
     * the genome Reads that overlap with the end of the genome, are currently
     * cut till the end of the genome and the overlap is deleted.
     *
     * @param input
     * @return
     * @throws java.io.IOException
     */
    private void inspectReads(SAMRecord input) throws IOException, CloneNotSupportedException {
        if(input.getReferenceName().equals("*")){
            output_realigned.addAlignment(input);
            return;
        }

        // We always keep un-mapped reads by ensuring they have a reference that is one of changedSequences.
        // We discard mapped reads if filter_sequence_ids or filter is set and they do not map to a reference
        // that is one of changedSequences.
        // If no filter is set and they do not map to a reference that is one of the changedSequences we keep
        // the reads.
        if(input.getReadUnmappedFlag()){
            Iterator<String> changedSequencesIterator = this.changedSequences.iterator ();
            if ( changedSequencesIterator.hasNext() && !this.changedSequences.contains(input.getReferenceName()) ) {
               String changedSequenceName = changedSequencesIterator.next ();
               input.setReferenceName (changedSequenceName);
               input.setReferenceIndex (this.modSamheader.getSequenceIndex(input.getReferenceName()));
               input.setAlignmentStart (1);
            }
            output_realigned.addAlignment(input);
            return;
	} else if( !this.changedSequences.contains(input.getReferenceName()) && (filter_sequence_ids || filtered) ){
          return;
        } else if( !this.changedSequences.contains(input.getReferenceName()) ) {
          output_realigned.addAlignment(input);
          return;
        }
        int startPos = input.getAlignmentStart();
        String seq = input.getReadString();

        // Useful variables
        int readLength = seq.length();

        //This is the most prominent, regular case!
        // Reads are within the non-modified area
        int elongatedReferenceLength = inputSam.getFileHeader().getSequence(input.getReferenceName()).getSequenceLength();
        int originalReferenceLength = elongatedReferenceLength - this.elongationfactor;

        if(startPos + readLength >= this.elongationfactor && startPos + readLength < originalReferenceLength) {
            stat_inregulargenome++;
            input.setReferenceIndex (this.modSamheader.getSequenceIndex(input.getReferenceName()));
            output_realigned.addAlignment(input);
        } else {
            processReads(input);
        }
    }

    private void processReads(SAMRecord read) throws CloneNotSupportedException {
        LinkedList<SAMRecord> tempSplit = new LinkedList<SAMRecord>();


            int startPos = read.getAlignmentStart();
            
            int elongatedReferenceLength = inputSam.getFileHeader().getSequence(read.getReferenceName()).getSequenceLength();
            int originalReferenceLength = elongatedReferenceLength - this.elongationfactor;

            //Then we're dealing with reads which could be potentially mapped twice, once at the beginning and once at the end of the mod genome, thus we should check the flags

            if(((startPos > originalReferenceLength) && (read.getReadLength() + startPos <= elongatedReferenceLength)) || startPos + read.getReadLength() <= this.elongationfactor) {
                stat_infactor++;

                String xa_tag = (String) read.getAttribute("XA");
                //If this holds true, we only have one duplicate read, which is truly no real duplicate -> just set MAPQ to 37 and get it back into alignment.
                if (xa_tag != null) {
                    if (xa_tag.endsWith(";")) {
                        read.setMappingQuality(37);
                        if(startPos > originalReferenceLength){
                            read.setAlignmentStart(returnMTPosition(xa_tag,read.getReferenceName()));
                        }
                        stat_inelongation++;
                        read.setAttribute("XA", ""); //Kill the attribute, the tag should be gone after determining if its a single hit!
                        read.setReadName(read.getReadName()+"C");
                        read.setReferenceIndex (this.modSamheader.getSequenceIndex(read.getReferenceName()));
                        read.setMateReferenceIndex(this.modSamheader.getSequenceIndex(read.getReferenceName()));
                        output_realigned.addAlignment(read);
                    }
                }
                else {
                    //then we need to set alignment start properly
                    read.setReadName(read.getReadName()+"C");
                    if(startPos > originalReferenceLength){
                        read.setAlignmentStart(startPos-originalReferenceLength);
                    } else {
                        //do nothing, its already mapped properly ;-)
                    }
                    read.setReferenceIndex (this.modSamheader.getSequenceIndex(read.getReferenceName()));
                    output_realigned.addAlignment(read);
                }
            }else {
            //Reads that lie in the intersection zone and thus need to be split into two parts
                splitRead(read, tempSplit);
                stat_overlapping++;
               // remove_reads.add(read);

            }


        flushMappedReads(tempSplit);
    }


    private void flushMappedReads(LinkedList<SAMRecord> list) {
        for(SAMRecord read : list){
            read.setReferenceIndex (this.modSamheader.getSequenceIndex(read.getReferenceName()));
            output_realigned.addAlignment(read);
        }
    }

    private void splitRead(SAMRecord read, LinkedList<SAMRecord> tempSplit) throws CloneNotSupportedException {
        SAMRecord copyBegin = (SAMRecord) read.clone();
        SAMRecord copyEnd = (SAMRecord) read.clone();

        int startPosition = read.getAlignmentStart();
        
        int elongatedReferenceLength = inputSam.getFileHeader().getSequence(read.getReferenceName()).getSequenceLength();
        int originalReferenceLength = elongatedReferenceLength - this.elongationfactor;

        int diff = 0;
        if(startPosition > originalReferenceLength){
             diff = elongatedReferenceLength - startPosition;
        } else {
             diff = originalReferenceLength - startPosition;

        }


        copyBegin.setReadName(read.getReadName()+"A");



        rt.trimRead(copyBegin, diff , true);
        if(startPosition < originalReferenceLength){
            copyBegin.setAlignmentStart(startPosition);
        } else {
            copyBegin.setAlignmentStart(startPosition-originalReferenceLength);
        }

        int splitDiff = Math.abs(read.getReadLength() - diff);

        copyEnd.setReadName(read.getReadName() + "B");
        rt.trimRead(copyEnd, read.getReadLength() - splitDiff +1, false);

        if(copyEnd.getCigar().isEmpty() || (copyEnd.getCigarLength() <=1 ) || copyEnd.getReadLength() == 0){

        } else {
            //This has to be set to 1, or else picard validateSamFile will complain about wrong start sites.
            copyEnd.setAlignmentStart(1);
            //add the splitted parts

            tempSplit.add(copyEnd);
        }
        
        //add the splitted parts

        if(copyBegin.getCigar().isEmpty() || (copyBegin.getCigarLength() <=1 ) || copyBegin.getReadLength() == 0){

        } else {
            tempSplit.add(copyBegin);

        }


    }

    /**
     * Method providing the CLI for the main method, also parsing input parameters etc.
     * @param args
     * @return
     */

    private static RealignSAMFile addCLIInterface(String[] args) {
        Options helpOptions = new Options();
        helpOptions.addOption("h", "help", false, "show this help page");
        Options options = new Options();
        options.addOption("h", "help", false, "show this help page");
        options.addOption(OptionBuilder.withLongOpt("input")
                .withArgName("INPUT")
                .withDescription("the input SAM/BAM File")
                .isRequired()
                .hasArg()
                .create("i"));
        options.addOption(OptionBuilder.withLongOpt("reference")
                .withArgName("REFERENCE")
                .withDescription("the unmodified reference genome")
                .isRequired()
                .hasArg()
                .create("r"));
        options.addOption(OptionBuilder.withLongOpt("elongation")
                .withArgName("ELONGATION")
                .withDescription("the elongation factor [INT]")
                .isRequired()
                .hasArg()
                .create("e"));
        options.addOption(OptionBuilder.withLongOpt("filterCircularReads")
                .isRequired(false)
                .withArgName("FILTER")
                .withDescription("filter the reads that don't map to a circular identifier (true/false), default false")
                .hasArg()
                .create("f"));
        options.addOption(OptionBuilder.withLongOpt("filterNonCircularSequences")
                .isRequired(false)
                .withArgName("FILTERSEQUENCEIDS")
                .withDescription("filter the sequence identifiers that are not circular identifiers (true/false), default false")
                .hasArg()
                .create("x"));
        HelpFormatter helpformatter = new HelpFormatter();
        CommandLineParser parser = new BasicParser();
        try {
            CommandLine cmd = parser.parse(helpOptions, args);
            if (cmd.hasOption('h')) {
                helpformatter.printHelp(CLASS_NAME +"v"+ VERSION, options);
                System.exit(0);
            }
        } catch (ParseException e1) {
        }
        String input = "";
        String reference = "";
        String tmpElongation = "";
        Integer elongation = 0;
        String filter = "";
	String filterSequenceIds = "";
        try {
            CommandLine cmd = parser.parse(options, args);
            if (cmd.hasOption('i')) {
                input = cmd.getOptionValue('i');
            }
            if (cmd.hasOption('e')) {
                tmpElongation = cmd.getOptionValue('e');
                try{
                    elongation = Integer.parseInt(tmpElongation);
                }catch(Exception e){
                    System.err.println("elongation not an Integer: " + tmpElongation);
                    System.exit(0);
                }
            }
            if(cmd.hasOption('r')){
                reference = cmd.getOptionValue('r');
            }
            if(cmd.hasOption('f')){
                filter = cmd.getOptionValue('f');
                if(filter.equals("true") || filter.equals("TRUE") || filter.equals("1")){
                    filtered = true;
                } else {
                    filtered = false;
                }
            }
            if (cmd.hasOption('x')) {
               filterSequenceIds = cmd.getOptionValue ('x');
               if ( filterSequenceIds.equals("true") || filter.equals("TRUE") || filter.equals("1")) {
                 filter_sequence_ids = true;
               } else {
                 filter_sequence_ids = false;
               }
            }
        } catch (ParseException e) {
            helpformatter.printHelp(CLASS_NAME, options);
            System.err.println(e.getMessage());
            System.exit(0);
        }

        RealignSAMFile rs = new RealignSAMFile(new File(input), new File(reference), elongation, filter);
        return rs;
    }


    private int returnMTPosition(String xa_tag, String toExtend){
        //e.g. XA:Z:chr5,-79948109,46M,1;chrMT,+374,46M,1;
        String[] splitMT = xa_tag.split(";");
        for(String s : splitMT){
            if(s.split(",")[0].equals(toExtend)){
                int index = Integer.parseInt(s.split(",")[1]);
                return Math.abs(index);
            }
        }
        return 1;

    }


}
