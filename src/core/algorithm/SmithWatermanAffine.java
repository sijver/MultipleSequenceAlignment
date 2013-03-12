package core.algorithm;

import core.Protein;
import core.SubstitutionMatrix;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class SmithWatermanAffine extends AlignmentAlgorithmAbstract implements AlignmentAlgorithm {

    private int[][] vMatrix;
    private int[][] wMatrix;
    private int extendingGap;

    public SmithWatermanAffine(Protein protein1, Protein protein2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        super(protein1, protein2, substitutionMatrix, gapPenalty);
    }

    public void computeSMatrix() {
        sMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        vMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        wMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];

        for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
            for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
                vMatrix[i][j] = Math.max(sMatrix[i - 1][j] + gapPenalty, vMatrix[i - 1][j] + extendingGap);
                wMatrix[i][j] = Math.max(sMatrix[i][j - 1] + gapPenalty, wMatrix[i][j - 1] + extendingGap);
                sMatrix[i][j] = Math.max(Math.max(vMatrix[i][j], wMatrix[i][j]), sMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j)));
                sMatrix[i][j] = Math.max(sMatrix[i][j], 0);
            }
        }
    }

    public void computeAlignments() {
        computeSMatrix();
        searchMatrixMaxValue(false);
        score = sMatrix[matrixMaxI][matrixMaxJ];

        List<String> alignment1List = new LinkedList<String>();
        List<String> alignment2List = new LinkedList<String>();

        int i = matrixMaxI;
        int j = matrixMaxJ;
        while (i > 0 && j > 0) {
            int score = sMatrix[i][j];
            int diagScore = sMatrix[i - 1][j - 1];
            int upScore = wMatrix[i][j];
            int leftScore = vMatrix[i][j];
            if (score == diagScore + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j))) {
                alignment1List.add(protein1.getProteinString().substring(i - 1, i));
                alignment2List.add(protein2.getProteinString().substring(j - 1, j));
                i--;
                j--;
            } else if(score == leftScore){
                alignment1List.add(protein1.getProteinString().substring(i - 1, i));
                alignment2List.add("-");
                i--;
            } else if (score == upScore) {
                alignment1List.add("-");
                alignment2List.add(protein2.getProteinString().substring(j - 1, j));
                j--;
            } else {
                break;
            }
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

    public void setExtendingGap(int extendingGap) {
        this.extendingGap = extendingGap;
    }
}
