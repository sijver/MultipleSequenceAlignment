package core.algorithm;

import core.Protein;
import core.SubstitutionMatrix;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class NeedlemanWunschAffine extends AlignmentAlgorithmAbstract implements AlignmentAlgorithm {

    private int[][] vMatrix;
    private int[][] wMatrix;
    private int extendingGap;
    private boolean isSemiGlobal = false;

    public NeedlemanWunschAffine(Protein protein1, Protein protein2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        super(protein1, protein2, substitutionMatrix, gapPenalty);
    }

    public void computeSMatrix() {
        sMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        vMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        wMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];

        //Matrix initialization
        if (!isSemiGlobal) {
            vMatrix[0][0] = -100000;
            wMatrix[0][0] = -100000;
            for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
                vMatrix[i][0] = gapPenalty + extendingGap * (i - 1);
                sMatrix[i][0] = -100000;
                wMatrix[i][0] = -100000;
            }
            for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
                wMatrix[0][j] = gapPenalty + extendingGap * (j - 1);
                sMatrix[0][j] = -100000;
                vMatrix[0][j] = -100000;
            }
        }

        for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
            for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
                vMatrix[i][j] = Math.max(sMatrix[i - 1][j] + gapPenalty, vMatrix[i - 1][j] + extendingGap);
                wMatrix[i][j] = Math.max(sMatrix[i][j - 1] + gapPenalty, wMatrix[i][j - 1] + extendingGap);
                sMatrix[i][j] = Math.max(Math.max(vMatrix[i - 1][j - 1], wMatrix[i - 1][j - 1]), sMatrix[i - 1][j - 1]) + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j));
            }
        }

        System.out.println("v");
        for(int[] i : vMatrix){
            System.out.println(Arrays.toString(i));
        }

        System.out.println("w");
        for(int[] i : wMatrix){
            System.out.println(Arrays.toString(i));
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
            int diagScore = sMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j));
            int upScore = wMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j));
            int leftScore = vMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j));
            if (score == diagScore) {
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

    public void setExtendingGap(int extendingGap) {
        this.extendingGap = extendingGap;
    }

    public void setSemiGlobal(boolean semiGlobal) {
        isSemiGlobal = semiGlobal;
    }
}
