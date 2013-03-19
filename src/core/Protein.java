package core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Protein {

    private String proteinString;
    private String id;

    public Protein(String id, String proteinString) {
        this(proteinString);
        this.id = id;
    }

    public Protein(String proteinString) {
        this.proteinString = proteinString.replaceAll("\\s", "");
    }

    public String getProteinString() {
        return this.proteinString;
    }

    public void setProteinString(String proteinString) {
        this.proteinString = proteinString;
    }

    public String getId() {
        return this.id;
    }

    public void setId(String id) {
        this.id = id;
    }

    @Override
    public String toString() {
        return "Dna [id=" + this.id + ", proteinString=" + this.proteinString
                + "]";
    }

    public static Protein shuffleProtein(Protein protein) {
        String proteinString = protein.getProteinString();
        StringBuilder shuffledProteinString = new StringBuilder(proteinString.length());
        List<Integer> acidPositions = new ArrayList<Integer>();
        for (int i = 0; i < proteinString.length(); i++) {
            acidPositions.add(i);
        }
        Collections.shuffle(acidPositions);
        for (Integer acidPosition : acidPositions) {
            shuffledProteinString.append(proteinString.substring(acidPosition, acidPosition + 1));
        }
        return new Protein(shuffledProteinString.toString());
    }

}