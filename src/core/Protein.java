package core;

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

}