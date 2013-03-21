package core.msa;

import core.Protein;
import core.SubstitutionMatrix;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class MultipleSequenceAlignment {

    private int openingGap;
    private int extendingGap;
    private SubstitutionMatrix substitutionMatrix;

    public MultipleSequenceAlignment(int openingGap, int extendingGap, SubstitutionMatrix substitutionMatrix) {
        this.openingGap = openingGap;
        this.extendingGap = extendingGap;
        this.substitutionMatrix = substitutionMatrix;
    }

    public List<Protein> makeMultipleSequenceAlignment(Cluster cluster) {
        List<Protein> proteinList1 = new LinkedList<Protein>();
        if (cluster.getClusterObjects().get(0).getClass() == Protein.class) {
            proteinList1.add((Protein) cluster.getClusterObjects().get(0));
        } else {
            proteinList1 = makeMultipleSequenceAlignment((Cluster) cluster.getClusterObjects().get(0));
        }
        List<Protein> proteinList2 = new LinkedList<Protein>();
        if (cluster.getClusterObjects().get(1).getClass() == Protein.class) {
            proteinList2.add((Protein) cluster.getClusterObjects().get(1));
        } else {
            proteinList2 = makeMultipleSequenceAlignment((Cluster) cluster.getClusterObjects().get(1));
        }
        return groupToGroupAlignment(proteinList1, proteinList2);
    }

    public List<Protein> groupToGroupAlignment(List<Protein> proteinList1, List<Protein> proteinList2) {
        int I = proteinList1.get(0).getProteinString().length();
        int J = proteinList2.get(0).getProteinString().length();
        int M = proteinList1.size();
        int N = proteinList2.size();
        int V = M * N * openingGap;
        double[][] memoryD1 = new double[I + 1][J + 1];
        double[][] memoryD2 = new double[I + 1][J + 1];
        double[][] memoryD3 = new double[I + 1][J + 1];
        for (int i = 1; i <= I; i++) {
            for (int j = 1; j <= J; j++) {
                memoryD1[i][j] = Math.min(memoryD3[i - 1][j] + V, memoryD1[i - 1][j]);
                int sumOfDistances = 0;
                for (int m = 0; m < M; m++) {
                    for (int n = 0; n < N; n++) {
                        sumOfDistances += aminoAcidDistance(proteinList1.get(m).getProteinString().substring(i - 1, i), "-");
                    }
                }
                memoryD1[i][j] += sumOfDistances;
                memoryD2[i][j] = Math.min(memoryD3[i][j - 1] + V, memoryD2[i][j - 1]);
                sumOfDistances = 0;
                for (int m = 0; m < M; m++) {
                    for (int n = 0; n < N; n++) {
                        sumOfDistances += aminoAcidDistance("-", proteinList2.get(n).getProteinString().substring(j - 1, j));
                    }
                }
                memoryD2[i][j] += sumOfDistances;
                memoryD3[i][j] = Math.min(Math.min(memoryD1[i - 1][j - 1], memoryD2[i - 1][j - 1]), memoryD3[i - 1][j - 1]);
                sumOfDistances = 0;
                for (int m = 0; m < M; m++) {
                    for (int n = 0; n < N; n++) {
                        sumOfDistances += aminoAcidDistance(proteinList1.get(m).getProteinString().substring(i - 1, i), proteinList2.get(n).getProteinString().substring(j - 1, j));
                    }
                }
                memoryD3[i][j] += sumOfDistances;
            }
        }
        List<Protein> listOfAlignedProteins = new LinkedList<Protein>();
        List<StringBuilder> listOfProteinAcids = new LinkedList<StringBuilder>();
        for (int i = 0; i < M + N; i++) {
            listOfProteinAcids.add(new StringBuilder());
        }

        int i = I;
        int j = J;
        while (i > 0 && j > 0) {
            if (memoryD1[i - 1][j - 1] <= memoryD2[i - 1][j - 1] && memoryD1[i - 1][j - 1] <= memoryD3[i - 1][j - 1]) {
                for (int k = 0; k < M + N; k++) {
                    if (k < M) {
                        listOfProteinAcids.get(k).append(proteinList1.get(k).getProteinString().substring(i - 1, i));
                    } else {
                        listOfProteinAcids.get(k).append("-");
                    }
                }
                i--;
            } else if (memoryD2[i - 1][j - 1] <= memoryD1[i - 1][j - 1] && memoryD2[i - 1][j - 1] <= memoryD3[i - 1][j - 1]) {
                for (int k = 0; k < M + N; k++) {
                    if (k < M) {
                        listOfProteinAcids.get(k).append("-");
                    } else {
                        listOfProteinAcids.get(k).append(proteinList2.get(k - M).getProteinString().substring(j - 1, j));
                    }
                }
                j--;
            } else {
                for (int k = 0; k < M + N; k++) {
                    if (k < M) {
                        listOfProteinAcids.get(k).append(proteinList1.get(k).getProteinString().substring(i - 1, i));
                    } else {
                        listOfProteinAcids.get(k).append(proteinList2.get(k - M).getProteinString().substring(j - 1, j));
                    }
                }
                i--;
                j--;
            }
        }
        while (i > 0) {
            for (int k = 0; k < M + N; k++) {
                if (k < M) {
                    listOfProteinAcids.get(k).append(proteinList1.get(k).getProteinString().substring(i - 1, i));
                } else {
                    listOfProteinAcids.get(k).append("-");
                }
            }
            i--;
        }
        while (j > 0) {
            for (int k = 0; k < M + N; k++) {
                if (k < M) {
                    listOfProteinAcids.get(k).append("-");
                } else {
                    listOfProteinAcids.get(k).append(proteinList2.get(k - M).getProteinString().substring(j - 1, j));
                }
            }
            j--;
        }

        for (int k = 0; k < listOfProteinAcids.size(); k++) {
            String proteinString = listOfProteinAcids.get(k).reverse().toString();
            String proteinId;
            if (k < M) {
                proteinId = proteinList1.get(k).getId();
            } else {
                proteinId = proteinList2.get(k - M).getId();
            }
            listOfAlignedProteins.add(new Protein(proteinId, proteinString));
        }
        return listOfAlignedProteins;
    }

    public int aminoAcidDistance(String acid1, String acid2) {
        if (acid1.equals("-") && acid2.equals("-")) {
            return 0;
        } else if (acid1.equals("-") || acid2.equals("-")) {
            return extendingGap;
        } else {
            return substitutionMatrix.getMatrixCell(acid1, acid2);
        }
    }

}
