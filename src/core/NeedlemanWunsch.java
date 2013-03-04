package core;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class NeedlemanWunsch {

    private Protein protein1;
    private Protein protein2;
    private int gapPenalty;
    private SubstitutionMatrix substitutionMatrix;
    private int[][] sMatrix;
    private int[][] vMatrix;
    private int[][] wMatrix;
    private StringBuilder alignment1;
    private StringBuilder alignment2;
    private double identity;
    private int identicalPositions;

    public NeedlemanWunsch(Protein protein1, Protein protein2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        this.protein1 = protein1;
        this.protein2 = protein2;
        this.substitutionMatrix = substitutionMatrix;
        this.gapPenalty = gapPenalty;
    }

    private int getAlignment(String s1, String s2) {
        return substitutionMatrix.getMatrixCell(s1, s2);
    }

    private void computeSMatrix() {
        sMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];

        //Matrix initialization
        for (int i = 0; i < protein1.getProteinString().length() + 1; i++) {
            sMatrix[i][0] = gapPenalty * i;
        }
        for (int j = 0; j < protein2.getProteinString().length() + 1; j++) {
            this.sMatrix[0][j] = gapPenalty * j;
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

    private void computeMatrices(int extendingGap) {
        sMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        vMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];
        wMatrix = new int[protein1.getProteinString().length() + 1][protein2.getProteinString().length() + 1];

        //Matrix initialization
        for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
            sMatrix[i][0] = gapPenalty + extendingGap * (i - 1);
            vMatrix[i][0] = gapPenalty + extendingGap * (i - 1);
            wMatrix[i][0] = gapPenalty + extendingGap * (i - 1);
        }
        for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
            sMatrix[0][j] = gapPenalty + extendingGap * (j - 1);
            vMatrix[0][j] = gapPenalty + extendingGap * (j - 1);
            wMatrix[0][j] = gapPenalty + extendingGap * (j - 1);
        }

        for (int i = 1; i < protein1.getProteinString().length() + 1; i++) {
            for (int j = 1; j < protein2.getProteinString().length() + 1; j++) {
                vMatrix[i][j] = Math.max(sMatrix[i - 1][j] + gapPenalty, vMatrix[i - 1][j] + extendingGap);
                wMatrix[i][j] = Math.max(sMatrix[i][j - 1] + gapPenalty, wMatrix[i][j - 1] + extendingGap);
                sMatrix[i][j] = Math.max(Math.max(vMatrix[i][j], wMatrix[i][j]), sMatrix[i - 1][j - 1] + getAlignment(protein1.getProteinString().substring(i - 1, i), protein2.getProteinString().substring(j - 1, j)));
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

    public String getAlignment1() {
        return alignment1.toString();
    }

    public String getAlignment2() {
        return alignment2.toString();
    }

    public void computeAlignments() {
        computeSMatrix();

        List<String> alignment1List = new LinkedList<String>();
        List<String> alignment2List = new LinkedList<String>();

        int i = protein1.getProteinString().length();
        int j = protein2.getProteinString().length();
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

    public void computeAlignmentsAffine(int extendingGap) {
        computeMatrices(extendingGap);

        List<String> alignment1List = new LinkedList<String>();
        List<String> alignment2List = new LinkedList<String>();

        int i = protein1.getProteinString().length();
        int j = protein2.getProteinString().length();
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

    private void computeIdentity(){
        identicalPositions = 0;
        for(int i = 0; i < alignment1.length(); i++){
            if(alignment1.charAt(i) == alignment2.charAt(i)){
                identicalPositions++;
            }
        }
        identity = (double)identicalPositions / alignment1.length();
    }

    public int getScore(){
        return sMatrix[sMatrix.length - 1][sMatrix[0].length - 1];
    }

    public int[][] getsMatrix() {
        return sMatrix;
    }

    public double getIdentity() {
        return identity;
    }

    public int getIdenticalPositions() {
        return identicalPositions;
    }
}
