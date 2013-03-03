package core.io;

import core.SubstitutionMatrix;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class SubstitutionMatrixReader {

    private String inputFile;
    private SubstitutionMatrix substitutionMatrix;

    public SubstitutionMatrixReader(String inputFile) {
        this.inputFile = inputFile;
    }

    public SubstitutionMatrix processFile() {

        BufferedReader fileReader = null;
        String line = null;
        String id = "";
        StringBuffer contentBuffer = new StringBuffer();

        boolean firstLine = true;
        int rowNum = 0;

        try {

            fileReader = new BufferedReader(new FileReader(inputFile));
            while ((line = fileReader.readLine()) != null) {
                    if(line.startsWith("#") || line.isEmpty()){
                        // comment line, skip it
                        continue;
                    }
                    if (firstLine) {
                        List<String> acids = Arrays.asList(line.trim().split("\\s+"));
                        substitutionMatrix = new SubstitutionMatrix(acids.size());
                        substitutionMatrix.setAminoacids(acids);
                        firstLine = false;
                    } else {
                        String[] tokens = line.trim().split("\\s+");
                        int columnNum = 0;
                        for(String token : tokens){
                            try {
                                substitutionMatrix.setMatrixCell(rowNum, columnNum, Integer.parseInt(token));
                                columnNum++;
                            } catch (NumberFormatException ignored){
                            }
                        }
                        rowNum++;
                    }

            }

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

        return substitutionMatrix;

    }

}
