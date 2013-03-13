package core.msa;

import core.Protein;
import core.SubstitutionMatrix;
import core.algorithm.NeedlemanWunschAffine;
import core.io.FastaReader;
import core.io.SubstitutionMatrixReader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class ProgressiveAlignment {

    private int[][] distanceMatrix;
    private int[][] randomMatrix;
    private double[][] resultingMatrix;
    private List<Cluster> clustersList;
    private List<Protein> proteins;
    private SubstitutionMatrix substitutionMatrix;
    private int closestClusterNum1;
    private int closestClusterNum2;

    public ProgressiveAlignment() {
    }

    public void readFastaFile(String fileName) {
        FastaReader fr = new FastaReader(fileName);
        proteins = fr.processFile();
        initClustersList();
    }

    public void readSubstitutionMatrix(String fileName) {
        SubstitutionMatrixReader smr = new SubstitutionMatrixReader(fileName);
        substitutionMatrix = smr.processFile();
    }

    public void calculateDistanceMatrix(int openingGap, int extendinGap) {
        distanceMatrix = new int[proteins.size()][proteins.size()];
        for (int i = 0; i < proteins.size(); i++) {
            for (int j = i; j < proteins.size(); j++) {
                //TODO algorithm choosing
                NeedlemanWunschAffine algorithm = new NeedlemanWunschAffine(proteins.get(i), proteins.get(j), substitutionMatrix, openingGap);
                algorithm.setExtendingGap(extendinGap);
                algorithm.computeAlignments();
                distanceMatrix[i][j] = algorithm.getScore();
                //TODO is needed?
                distanceMatrix[j][i] = algorithm.getScore();
            }
        }
    }

    public void calculateRandomMatrix(int openingGap, int extendinGap) {
        randomMatrix = new int[proteins.size()][proteins.size()];
        for (int i = 1; i < proteins.size(); i++) {
            for (int j = i + 1; j < proteins.size(); j++) {
                //TODO algorithm choosing
                NeedlemanWunschAffine algorithm = new NeedlemanWunschAffine(shuffleProtein(proteins.get(i)), shuffleProtein(proteins.get(j)), substitutionMatrix, openingGap);
                algorithm.setExtendingGap(extendinGap);
                algorithm.computeAlignments();
                randomMatrix[i][j] = algorithm.getScore();
                //TODO is needed?
                randomMatrix[j][i] = algorithm.getScore();
            }
        }
    }

    public void calculateResultingMatrix() {
        resultingMatrix = new double[proteins.size()][proteins.size()];
        for (int i = 1; i < proteins.size(); i++) {
            for (int j = i + 1; j < proteins.size(); j++) {
                resultingMatrix[i][j] = -Math.log((distanceMatrix[i][j] - randomMatrix[i][j]) / ((randomMatrix[i][i] + randomMatrix[j][j]) / 2 - randomMatrix[i][j]));
            }
        }
    }

    public void initClustersList() {
        clustersList = new LinkedList<Cluster>();
        for (int i = 0; i < proteins.size(); i++) {
            Cluster newCLuster = new Cluster();
            newCLuster.addObjectToCluster(i);
            clustersList.add(newCLuster);
        }
    }

    public void findSmallestDistance() {
        int minDistance = distanceMatrix[0][1];
        closestClusterNum1 = 0;
        closestClusterNum2 = 1;
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i; j < distanceMatrix[0].length; j++) {
                if (distanceMatrix[i][j] <= minDistance) {
                    minDistance = distanceMatrix[i][j];
                    closestClusterNum1 = i;
                    closestClusterNum2 = j;
                }
            }
        }
    }

    public void make() {
        findSmallestDistance();
        Cluster cluster1 = clustersList.get(closestClusterNum1);
        Cluster cluster2 = clustersList.get(closestClusterNum2);
        Cluster cluster3 = new Cluster();
        for (int i = 0; i < clustersList.size(); i++) {
            if (i != closestClusterNum1 && i != closestClusterNum2) {
                cluster3.addObjectToCluster(clustersList.get(i));
            }
        }

    }

    public Protein shuffleProtein(Protein protein) {
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
