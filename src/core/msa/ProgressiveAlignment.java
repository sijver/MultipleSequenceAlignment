package core.msa;

import core.Protein;
import core.SubstitutionMatrix;
import core.algorithm.*;
import core.io.FastaReader;
import core.io.SubstitutionMatrixReader;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class ProgressiveAlignment {

    private List<Protein> proteins; //List of proteins to align
    private AlignmentAlgorithm alignmentAlgorithm;  //Algorithm used for two sequences alignment
    private AlgorithmType algorithmType;    //Type of algorithm used for alignment
    private SubstitutionMatrix substitutionMatrix;  //Substitution matrix used for sequence alignment
    private int openingGap; //Opening gap used for sequence alignment
    private int extendingGap;   //Extending gap used for sequence alignment
    private int[][] alignmentMatrix;
    private int[][] randomMatrix;
    private double[][] distanceMatrix;
    private List<Cluster> clustersList;

    public ProgressiveAlignment(List<Protein> proteinsList, AlgorithmType algorithmType, SubstitutionMatrix substitutionMatrix, int openingGap, int extendingGap) {
        this.proteins = proteinsList;
        this.algorithmType = algorithmType;
        this.substitutionMatrix = substitutionMatrix;
        this.openingGap = openingGap;
        this.extendingGap = extendingGap;
        initClustersList();
        calculateAlignmentMatrix();
        calculateRandomMatrix();
        calculateDistanceMatrix();
    }

    private void initClustersList() {
        clustersList = new LinkedList<Cluster>();
        for (int i = 0; i < proteins.size(); i++) {
            Cluster newCluster = new Cluster();
            newCluster.addObjectToCluster(i);
            clustersList.add(newCluster);
        }
    }

    private int proteinsAlignment(Protein protein1, Protein protein2){
        initAlgorithm(protein1, protein2);
        alignmentAlgorithm.computeAlignments();
        return alignmentAlgorithm.getScore();
    }

    private void initAlgorithm(Protein protein1, Protein protein2){
        switch (algorithmType){
            case NEEDLEMAN_WUNSCH_LINEAR:
                alignmentAlgorithm = new NeedlemanWunschLinear(protein1, protein2, substitutionMatrix, openingGap);
                break;
            case NEEDLEMAN_WUNSCH_AFFINE:
                alignmentAlgorithm = new NeedlemanWunschAffine(protein1, protein2, substitutionMatrix, openingGap);
                ((NeedlemanWunschAffine)alignmentAlgorithm).setExtendingGap(extendingGap);
                break;
            case SMITH_WATERMAN_LINEAR:
                alignmentAlgorithm = new SmithWatermanLinear(protein1, protein2, substitutionMatrix, openingGap);
                break;
            case SMITH_WATERMAN_AFFINE:
                alignmentAlgorithm = new SmithWatermanAffine(protein1, protein2, substitutionMatrix, openingGap);
                ((SmithWatermanAffine)alignmentAlgorithm).setExtendingGap(extendingGap);
        }
    }

    // Calculates alignment between all the the clusters
    private void calculateAlignmentMatrix() {
        alignmentMatrix = new int[proteins.size()][proteins.size()];
        for (int i = 0; i < proteins.size(); i++) {
            for (int j = i; j < proteins.size(); j++) {
                alignmentMatrix[i][j] = proteinsAlignment(proteins.get(i), proteins.get(j));
            }
        }
    }

    private void calculateRandomMatrix() {
        randomMatrix = new int[proteins.size()][proteins.size()];
        for (int i = 0; i < proteins.size(); i++) {
            for (int j = i + 1; j < proteins.size(); j++) {
                randomMatrix[i][j] = proteinsAlignment(Protein.shuffleProtein(proteins.get(i)), Protein.shuffleProtein(proteins.get(j)));
            }
        }
    }

    private void calculateDistanceMatrix() {
        distanceMatrix = new double[proteins.size()][proteins.size()];
        for (int i = 0; i < proteins.size(); i++) {
            for (int j = i + 1; j < proteins.size(); j++) {
                distanceMatrix[i][j] = -Math.log(((double)alignmentMatrix[i][j] - (double)randomMatrix[i][j]) / (((double)alignmentMatrix[i][i] + (double)alignmentMatrix[j][j]) / (double)2 - (double)randomMatrix[i][j]));
                distanceMatrix[j][i] = distanceMatrix[i][j];
            }
        }
    }

    private double clusterClusterDistance(Cluster cluster1, Cluster cluster2){
        List<Integer> clusterProteins1 = Cluster.getClusterProteins(cluster1);
        List<Integer> clusterProteins2 = Cluster.getClusterProteins(cluster2);
        double clustersDistance = Double.MAX_VALUE;
        for(int proteinNumber1 : clusterProteins1){
            for(int proteinNumber2 : clusterProteins2){
                double proteinDistance = distanceMatrix[proteinNumber1][proteinNumber2];
                if(proteinDistance < clustersDistance){
                    clustersDistance = proteinDistance;
                }
            }
        }
        return clustersDistance;
    }

    public void constructEvolutionaryTree() {
        int clusterMinNum1 = 0;
        int clusterMinNum2 = 1;
        double clustersMinDistance;
        double clustersDistance;
        while(clustersList.size() > 1){
            clustersMinDistance = Double.MAX_VALUE;
            for(int i = 0; i < clustersList.size(); i++){
                for(int j = i + 1; j < clustersList.size(); j++){
                    clustersDistance = clusterClusterDistance(clustersList.get(i), clustersList.get(j));
                    System.out.println("cld "+i+ "  "+j+"   "+clustersDistance);
                    if(clustersDistance < clustersMinDistance){
                        clusterMinNum1 = i;
                        System.out.println(clusterMinNum1+"aaaaaaaa");
                        clusterMinNum2 = j;
                        clustersMinDistance = clustersDistance;
                    }
                }
            }
            System.out.println(clustersList.size()+"size");
            System.out.println(clusterMinNum1);
            System.out.println(clusterMinNum2);
            Cluster newCluster = new Cluster();
            newCluster.addObjectToCluster(clustersList.get(clusterMinNum1));
            newCluster.addObjectToCluster(clustersList.get(clusterMinNum2));
            clustersList.remove(clustersList.get(clusterMinNum1));
            clustersList.remove(clustersList.get(clusterMinNum2 - 1));
            clustersList.add(newCluster);
        }
    }

    public List<Cluster> getClustersList() {
        return clustersList;
    }

    public static void main(String[] args) {
        FastaReader fr = new FastaReader("fasta.txt");
        List<Protein> proteins= fr.processFile();

        SubstitutionMatrixReader smr = new SubstitutionMatrixReader("pam250.bla");
        SubstitutionMatrix sm = smr.processFile();

        ProgressiveAlignment pa = new ProgressiveAlignment(proteins, AlgorithmType.NEEDLEMAN_WUNSCH_LINEAR, sm, -8, -1);
        pa.constructEvolutionaryTree();

        System.out.println(Cluster.getClusterTree(pa.getClustersList().get(0), 0));

    }

}
