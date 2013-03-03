package core.io;

import core.Protein;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */

public class FastaReader {

    private String inputFile;
    private List<Protein> proteinList;

    public FastaReader(String inputFile) {
        this.inputFile = inputFile;
        this.proteinList = new ArrayList<Protein>();
    }

    public List<Protein> processFile() {

        BufferedReader fileReader = null;
        String line = null;
        String id = "";
        StringBuffer contentBuffer = new StringBuffer();

        try {

            fileReader = new BufferedReader(new FileReader(inputFile));
            do {

                line = fileReader.readLine();
                if (line != null) {

                    line = line.trim();
                    if (line.isEmpty()) {
                        continue;
                    }
                    char firstChar = line.charAt(0);
                    if (firstChar == '>') {

                        // save the previous sequence read
                        addToProteinList(id, contentBuffer);

                        // now can get the new id > ..
                        id = line.substring(1).trim();

                        // start a new content buffer
                        contentBuffer = new StringBuffer();

                    } else if (firstChar == ';') {

                        // comment line, skip it

                    } else {

                        // carry on reading sequence content
                        contentBuffer.append(line.trim());

                    }

                } else {

                    // save the final sequence content
                    addToProteinList(id, contentBuffer);

                }

            } while (line != null);

        } catch (FileNotFoundException e) {
            System.out.println("File " + inputFile + " is not found");
            System.exit(1);
        } catch (IOException e) {
            System.out.println("An IO error has occured: " + e.getMessage());
            System.exit(1);
        } finally {
            if (fileReader != null) {
                try {
                    fileReader.close();
                } catch (IOException e) {
                    // do nothing
                }
            }
        }

        return proteinList;

    }

    private void addToProteinList(String id, StringBuffer sb) {
        if (sb.length() != 0) {
            Protein s = null;
            String content = sb.toString();
            s = new Protein(id, content);
            this.proteinList.add(s);
        }
    }

}