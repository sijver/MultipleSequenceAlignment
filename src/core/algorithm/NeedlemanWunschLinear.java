package core.algorithm;

import core.Protein;
import core.SubstitutionMatrix;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class NeedlemanWunschLinear extends AlignmentAlgorithmAbstract implements AlignmentAlgorithm {

    private boolean isSemiGlobal = false;

    public NeedlemanWunschLinear(Protein protein1, Protein protein2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        super(protein1, protein2, substitutionMatrix, gapPenalty);
    }

    public void computeSMatrix() {
        sMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];

        //Matrix initialization
        if (!isSemiGlobal) {
            for (int i = 0; i < protein1.getProteinString().length() + 1; i++) {
                sMatrix[i][0] = gapPenalty * i;
            }
            for (int j = 0; j < protein2.getProteinString().length() + 1; j++) {
                this.sMatrix[0][j] = gapPenalty * j;
            }
        }

        int match;
        int delete;
        int insert;
        for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
            for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
                match = sMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j));
                delete = sMatrix[i - 1][j] + gapPenalty;
                insert = sMatrix[i][j - 1] + gapPenalty;
                sMatrix[i][j] = Math.max(match, Math.max(delete, insert));
            }
        }
    }

    public void computeAlignments() {
        computeSMatrix();

        if (!isSemiGlobal) {
            score = sMatrix[sMatrix.length - 1][sMatrix[0].length - 1];
        } else {
            searchMatrixMaxValue(true);
            score = sMatrix[matrixMaxI][matrixMaxJ];
        }

        List<String> alignment1List = new LinkedList<String>();
        List<String> alignment2List = new LinkedList<String>();

        int i = protein1.getProteinString().length();
        int j = protein2.getProteinString().length();
        if(isSemiGlobal){
            while(i > matrixMaxI){
                alignment1List.add(protein1.getProteinString().substring(i - 1, i));
                alignment2List.add("-");
                i--;
            }
            while(j > matrixMaxJ){
                alignment1List.add("-");
                alignment2List.add(protein2.getProteinString().substring(j - 1, j));
                j--;
            }
        }

        while (i > 0 && j > 0) {
            int score = sMatrix[i][j];
            int diagScore = sMatrix[i - 1][j - 1];
            int upScore = sMatrix[i][j - 1];
            int leftScore = sMatrix[i - 1][j];
            if (score == diagScore + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j))) {
                alignment1List.add(protein1.getProteinString().substring(i - 1, i));
                alignment2List.add(protein2.getProteinString().substring(j - 1, j));
                i--;
                j--;
            } else if (score == leftScore + gapPenalty) {
                alignment1List.add(protein1.getProteinString().substring(i - 1, i));
                alignment2List.add("-");
                i--;
            } else {
                alignment1List.add("-");
                alignment2List.add(protein2.getProteinString().substring(j - 1, j));
                j--;
            }
        }
        while (i > 0) {
            alignment1List.add(protein1.getProteinString().substring(i - 1, i));
            alignment2List.add("-");
            i--;
        }
        while (j > 0) {
            alignment1List.add("-");
            alignment2List.add(protein2.getProteinString().substring(j - 1, j));
            j--;
        }

        // format alignments
        alignment1 = new StringBuilder(alignment1List.size());
        alignment2 = new StringBuilder(alignment2List.size());
        for (int k = alignment1List.size() - 1; k >= 0; k--) {
            alignment1.append(alignment1List.get(k));
        }
        for (int k = alignment2List.size() - 1; k >= 0; k--) {
            alignment2.append(alignment2List.get(k));
        }

        computeIdentity();

    }

    public void setSemiGlobal(boolean semiGlobal) {
        isSemiGlobal = semiGlobal;
    }
}
