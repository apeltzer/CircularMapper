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

import htsjdk.samtools.SAMRecord;

/**
 * Created by peltzer on 28/10/13.
 */
public class ReadTrimmer {


    public ReadTrimmer() {
    }

    public void trimRead(SAMRecord read, int trimLength, boolean fromStart) {


        String sequence = read.getReadString();
        String quality = read.getBaseQualityString();
        String cigar = read.getCigarString();

        String new_seq = "";
        String new_quali = "";
        String new_ziggi = ""; //Lulz

        // First String array = Integers between IDs.
        // Second String array = IDs
        String[] cigSplit = cigar.split("\\D");
        String[] intSplit = cigar.split("\\d+");

        if(fromStart){ //Trim Beginning

            trimLength++;
            int pos = 0;

            for(int i = 0; i< cigSplit.length; i++) {

                int tmp_int =  Integer.valueOf(cigSplit[i]);

                if(intSplit[i+1].equals("D")) {

                    trimLength = trimLength - tmp_int;

                    continue;

                } else if(intSplit[i+1].equals("I") ) {

                    new_seq += sequence.substring(pos,pos+tmp_int);
                    new_quali += quality.substring(pos,pos+tmp_int);

                    continue;

                } else if(intSplit[i+1].equals("S")){
                    new_seq += sequence.substring(pos,pos+tmp_int);
                    new_quali += quality.substring(pos,pos+tmp_int);
                    pos = pos+tmp_int;
                } else if(tmp_int < trimLength ){

                    trimLength = trimLength -  tmp_int;

                    new_seq += sequence.substring(pos, pos + tmp_int);
                    new_quali += quality.substring(pos, pos + tmp_int);

                    pos = pos + tmp_int;

                } else if (tmp_int >= trimLength ){

                    new_seq += sequence.substring(pos, pos + trimLength);
                    new_quali += quality.substring(pos, pos + trimLength);

                    cigSplit[i] = String.valueOf(trimLength);

                    for(int j = i +1; j< cigSplit.length; j++){
                        cigSplit[j] = "0";
                    }

                    break;
                }

            }


        } else { //Trim End

            int pos = sequence.length();
            trimLength = pos - trimLength;



            for(int i = cigSplit.length -1; i >= 0; i--) {

                int tmp_int =  Integer.valueOf(cigSplit[i]);

                if(intSplit[i+1].equals("D")){

                    trimLength = trimLength - tmp_int +1;
                    continue;

                }else if(intSplit[i+1].equals("I")) {

                    new_seq += sequence.substring(pos, pos + tmp_int);
                    new_quali += quality.substring(pos, pos + tmp_int);

                    continue;
                } else if(intSplit[i+1].equals("S")){
                    new_seq += sequence.substring(pos,pos+tmp_int);
                    new_quali += quality.substring(pos,pos+tmp_int);
                } else if(tmp_int < trimLength ){

                    trimLength = trimLength -  tmp_int;

                    new_seq = sequence.substring(pos - tmp_int,pos) + new_seq;
                    new_quali = quality.substring(pos - tmp_int,pos) + new_quali;

                    pos = pos - tmp_int;


                } else if (tmp_int > trimLength ){

                    new_seq = sequence.substring(pos - trimLength,pos) + new_seq;
                    new_quali = quality.substring(pos - trimLength,pos) + new_quali;

                    pos = pos - tmp_int;


                    cigSplit[i] = String.valueOf(trimLength);

                    for(int j = i -1; j>= 0; j--){
                        cigSplit[j] = "0";
                    }

                    break;
                }

            }



        }

        // Reconstruct CIGAR String
        for (int i = 0; i < cigSplit.length; i++) {
            if (Integer.parseInt(cigSplit[i]) >= 0) { // then we keep the CIGAR

                if (Integer.parseInt(cigSplit[i]) == 0) {
                    continue;
                }


                if(fromStart){
                    // String part
                    if( !(new_ziggi.equals("") &&  intSplit[i + 1].equals("D"))) {
                        new_ziggi += Integer.parseInt(cigSplit[i]) + intSplit[i + 1];
                    }
                } else {
                    // String part
                    if( !(new_ziggi.equals("") &&  intSplit[i + 1].equals("D"))) {
                        new_ziggi += Integer.parseInt(cigSplit[i]) + intSplit[i + 1];
                    }
                }

            } else {
                break;
            }
        }

        int diff = quality.length() - new_quali.length();
        if((diff > 0) & (diff < sequence.length() )){ //to get rid of reads with only masked part...
            if(!fromStart){
                new_ziggi = String.valueOf(diff) + "S" + new_ziggi;
            } else {
                new_ziggi = new_ziggi + String.valueOf(diff) + "S";
            }
            read.setCigarString(new_ziggi);
        } else {
            read.setCigarString(new_ziggi);
            read.setReadString("");
            read.setBaseQualityString("");
        }
    }
}
