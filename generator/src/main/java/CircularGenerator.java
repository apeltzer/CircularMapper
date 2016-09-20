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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Created with IntelliJ IDEA.
 * User: alex
 * Date: 27.08.13
 * Time: 11:13
 * To change this template use File | Settings | File Templates.
 */
public class CircularGenerator {
	private static final String CLASS_NAME = "CircularGenerator";
	private static final String VERSION = "1.0";
    private int elongationfactor = 0;
    private HashSet<String> keys_to_treat_circular = new HashSet<String>();

    public CircularGenerator(int elongationfactor) {
        this.elongationfactor = elongationfactor;
    }

    @SuppressWarnings("static-access")
	public static void main(String[] args) throws Exception {
    	Options helpOptions = new Options();
		helpOptions.addOption("h", "help", false, "show this help page");
		Options options = new Options();
		options.addOption("h", "help", false, "show this help page");
		options.addOption(OptionBuilder.withLongOpt("input")
				.withArgName("INPUT")
				.withDescription("the input FastA File")
				.isRequired()
				.hasArg()
				.create("i"));
		options.addOption(OptionBuilder.withLongOpt("elongation")
				.withArgName("ELONGATION")
				.withDescription("the elongation factor [INT]")
				.isRequired()
				.hasArg()
				.create("e"));
		options.addOption(OptionBuilder.withLongOpt("seq")
				.withArgName("SEQ")
				.withDescription("the names of the sequences that should to be elongated")
				.isRequired()
				.hasArg()
				.hasOptionalArgs()
				.hasArg()
				.create("s"));
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
		String tmpElongation = "";
		Integer elongation = 0;
		String[] names = new String[0];
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
			if (cmd.hasOption('s')) {
				names = cmd.getOptionValues('s');
			}
		} catch (ParseException e) {
			helpformatter.printHelp(CLASS_NAME, options);
			System.err.println(e.getMessage());
			System.exit(0);
		}
		
		CircularGenerator cg = new CircularGenerator(elongation);
		File f = new File(input);
		for(String s: names){
			cg.keys_to_treat_circular.add(s);
		}
		cg.extendFastA(f);



    }
    
    private void extendFastA(File inputFasta){
    	// create the new files
    	String fileExtension = com.google.common.io.Files.getFileExtension(inputFasta.getAbsolutePath());
        String fileName = com.google.common.io.Files.getNameWithoutExtension(inputFasta.getAbsolutePath());
        String filePath = inputFasta.getAbsolutePath().substring(0, inputFasta.getAbsolutePath().lastIndexOf(File.separator));
        File output = new File(filePath+"/"+fileName+"_"+this.elongationfactor+"."+fileExtension);
        File outputChanged = new File(inputFasta.getAbsoluteFile()+"_"+this.elongationfactor+"_elongated");
        PrintWriter out;
        // read the sequences
		try {
			out = new PrintWriter(new BufferedWriter(new FileWriter(output, false)));
			PrintWriter outChanged = new PrintWriter(new BufferedWriter(new FileWriter(outputChanged, false)));
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(inputFasta));
			String currLine = "";
			String header = "";
			Boolean firstEntry = true;
			StringBuffer sequence = new StringBuffer();
			// read the fasta file entry by entry
			while((currLine = br.readLine()) != null){
				// start of new sequence
				if(currLine.startsWith(">")){
					if(!firstEntry){
						extendSequence(header, sequence, out, outChanged);
					}
					// reset sequence
					firstEntry = false;
					header = currLine.substring(1);
					sequence.setLength(0);
				}else{
					sequence.append(currLine);
				}
			}
			// modify last entry
			extendSequence(header, sequence, out, outChanged);
			// close the writers
			out.flush();
			outChanged.flush();
			out.close();
			outChanged.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }


    private void extendSequence(String header, StringBuffer sequence, PrintWriter outFastA, PrintWriter outChanged) {
    	String currEntry = header.split(" ")[0]; //just use the forward part of the entry, never something else
    	//Modify the data
		boolean regEx = false;
		for(String names: keys_to_treat_circular){
			if(currEntry.contains(names)){
				regEx = true;
			}
		}
    	if(regEx){
    		if(sequence.length() < elongationfactor){
    			System.out.println(sequence.length());
                System.err.println("You cannot extend your sequence with a value longer than your actual sequence.");
                System.exit(1);
            } else {
				String extension = sequence.substring(0, this.elongationfactor);
				sequence.append(extension);
				outChanged.println(header.split(" ")[0]);
			}
    	}
    	// write the sequence to the new file
    	outFastA.println(">"+header);
    	int index = 0;
    	while(index <= sequence.length()){
    		outFastA.println(sequence.substring(index, Math.min(index+80, sequence.length())));
    		index += 80;
    	}
    	outFastA.flush();
	}


}
