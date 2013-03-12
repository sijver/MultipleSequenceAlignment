package core.algorithm;

import core.Protein;
import core.SubstitutionMatrix;

/**
 * Created with IntelliJ IDEA.
 */
public abstract class AlignmentAlgorithmAbstract implements AlignmentAlgorithm {

    Protein protein1;
    Protein protein2;
    int gapPenalty;
    SubstitutionMatrix substitutionMatrix;
    int[][] sMatrix;
    StringBuilder alignment1;
    StringBuilder alignment2;
    double identity;
    int identicalPositions;
    int score;

    int matrixMaxI = 0;
    int matrixMaxJ = 0;

    public AlignmentAlgorithmAbstract(Protein protein1, Protein protein2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        this.protein1 = protein1;
        this.protein2 = protein2;
        this.substitutionMatrix = substitutionMatrix;
        this.gapPenalty = gapPenalty;
    }

    public void computeIdentity() {
        identicalPositions = 0;
        for (int i = 0; i < alignment1.length(); i++) {
            if (alignment1.charAt(i) == alignment2.charAt(i)) {
                identicalPositions++;
            }
        }
        identity = (double) identicalPositions / alignment1.length();
    }

    public void searchMatrixMaxValue(boolean searchOverLastColumnAndRoOnly) {
        int maxValue = Integer.MIN_VALUE;
        if (!searchOverLastColumnAndRoOnly) {
            for (int i = 0; i < sMatrix.length; i++) {
                for (int j = 0; j < sMatrix[0].length; j++) {
                    if (sMatrix[i][j] > maxValue) {
                        maxValue = sMatrix[i][j];
                        matrixMaxI = i;
                        matrixMaxJ = j;
                    }
                }
            }
        } else {
            int matrixHeight = sMatrix.length;
            int matrixWidth = sMatrix[0].length;
            for (int i = 0; i < matrixHeight; i++) {
                if (sMatrix[i][matrixWidth - 1] > maxValue) {
                    maxValue = sMatrix[i][matrixWidth - 1];
                    matrixMaxI = i;
                    matrixMaxJ = matrixWidth - 1;
                }
            }
            for (int j = 0; j < matrixWidth; j++) {
                if (sMatrix[matrixHeight - 1][j] > maxValue) {
                    maxValue = sMatrix[matrixHeight - 1][j];
                    matrixMaxI = matrixHeight - 1;
                    matrixMaxJ = j;
                }
            }
        }
    }

    public int getAlignment(String s1, String s2) {
        return substitutionMatrix.getMatrixCell(s1, s2);
    }

    public int getScore() {
        return score;
    }

    public int[][] getSMatrix() {
        return sMatrix;
    }

    public double getIdentity() {
        return identity;
    }

    public int getIdenticalPositions() {
        return identicalPositions;
    }

    public String getAlignment1() {
        return alignment1.toString();
    }

    public String getAlignment2() {
        return alignment2.toString();
    }

}
