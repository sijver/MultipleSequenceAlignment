package core;

import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class SubstitutionMatrix {

    private int[][] substitutionMatrix;

    private List<String> aminoacids;

    public SubstitutionMatrix(int matrixSize) {
        substitutionMatrix = new int[matrixSize][matrixSize];
        aminoacids = new ArrayList<String>(matrixSize);
    }

    public void setMatrixCell(int i, int j, int value){
        substitutionMatrix[i][j] = value;
        substitutionMatrix[j][i] = value;
    }

    public int getMatrixCell(int i, int j){
        return substitutionMatrix[i][j];
    }

    public int getMatrixCell(String s1, String s2){
        return getMatrixCell(aminoacids.indexOf(s1), aminoacids.indexOf(s2));
    }

    public List<String> getAminoacids() {
        return aminoacids;
    }

    public void setAminoacids(List<String> aminoacids) {
        this.aminoacids = aminoacids;
    }

    public int[][] getSubstitutionMatrix() {
        return substitutionMatrix;
    }

}
