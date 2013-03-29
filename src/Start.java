import core.Protein;
import core.SubstitutionMatrix;
import core.algorithm.*;
import core.io.FastaReader;
import core.io.SubstitutionMatrixReader;
import core.msa.Cluster;
import core.msa.MultipleSequenceAlignment;
import core.msa.ProgressiveAlignment;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Start {

    public static void main(String[] args) {
        List<Protein> proteins;
        FastaReader fr = new FastaReader("fasta.txt");
        SubstitutionMatrix sm;
        SubstitutionMatrixReader smr = new SubstitutionMatrixReader("blosum62.bla");
        boolean isSemiGlobal = false;
        int openingGap = 10;
        int extendingGap = 1;
        boolean isMSA = true;
        AlignmentAlgorithm algorithm = null;
        AlgorithmType chosenAlgorithm = AlgorithmType.NEEDLEMAN_WUNSCH_LINEAR;


        if (args.length % 2 != 0) {
            wrongCommandExit();
        } else {
            for (int i = 0; i < args.length; i++) {
                if (i % 2 == 0) {
                    if (args[i].equals("-f")) {
                        fr = new FastaReader(args[i + 1]);
                    } else if (args[i].equals("-sm")) {
                        smr = new SubstitutionMatrixReader(args[i + 1]);
                    } else if (args[i].equals("-sg")) {
                        try {
                            if (Integer.parseInt(args[i + 1]) == 0) {
                                isSemiGlobal = false;
                            } else if (Integer.parseInt(args[i + 1]) == 1) {
                                isSemiGlobal = true;
                            } else {
                                wrongCommandExit();
                            }
                        } catch (Exception e) {
                            wrongCommandExit();
                        }
                    } else if (args[i].equals("-og")) {
                        try {
                            openingGap = Math.abs(Integer.parseInt(args[i + 1]));
                        } catch (Exception e) {
                            wrongCommandExit();
                        }
                    } else if (args[i].equals("-eg")) {
                        try {
                            extendingGap = Math.abs(Integer.parseInt(args[i + 1]));
                        } catch (Exception e) {
                            wrongCommandExit();
                        }
                    } else if (args[i].equals("-al")) {
                        if (args[i + 1].equals("msa")) {
                            isMSA = true;
                        } else if (args[i + 1].equals("tsa")) {
                            isMSA = false;
                        } else {
                            wrongCommandExit();
                        }
                    } else if (args[i].equals("-tsa")) {
                        if (args[i + 1].equals("nwl")) {
                            chosenAlgorithm = AlgorithmType.NEEDLEMAN_WUNSCH_LINEAR;
                        } else if (args[i + 1].equals("nwa")) {
                            chosenAlgorithm = AlgorithmType.NEEDLEMAN_WUNSCH_AFFINE;
                        } else if (args[i + 1].equals("swl")) {
                            chosenAlgorithm = AlgorithmType.SMITH_WATERMAN_LINEAR;
                        } else if (args[i + 1].equals("swa")) {
                            chosenAlgorithm = AlgorithmType.SMITH_WATERMAN_AFFINE;
                        } else {
                            wrongCommandExit();
                        }
                    } else {
                        wrongCommandExit();
                    }
                }
            }
        }
        proteins = fr.processFile();
        sm = smr.processFile();

        if (proteins.size() < 2) {
            System.out.println("\nFasta file contains less than 2 sequences\n");
            System.exit(1);
        }

        if (isMSA) {
            ProgressiveAlignment progressiveAlignment = new ProgressiveAlignment(proteins, AlgorithmType.NEEDLEMAN_WUNSCH_LINEAR, sm, openingGap * -1, extendingGap * -1);
            progressiveAlignment.constructEvolutionaryTree();
            System.out.println("\nGuide tree: ");
            System.out.println(Cluster.getClusterTree(progressiveAlignment.getRootCluster(), 0));

            System.out.println("\nAlignment: ");
            MultipleSequenceAlignment msa = new MultipleSequenceAlignment(openingGap, extendingGap, sm);
            List<Protein> proteins1 = msa.makeMultipleSequenceAlignment(progressiveAlignment.getRootCluster());
            for (Protein prot : proteins1) {
                System.out.println(prot.getProteinString() + " (" + prot.getId() + ")");
            }
            System.out.println();
        } else {
            if (chosenAlgorithm == AlgorithmType.NEEDLEMAN_WUNSCH_LINEAR) {
                algorithm = new NeedlemanWunschLinear(proteins.get(0), proteins.get(1), sm, openingGap * -1);
                ((NeedlemanWunschLinear) algorithm).setSemiGlobal(isSemiGlobal);
            } else if (chosenAlgorithm == AlgorithmType.NEEDLEMAN_WUNSCH_AFFINE) {
                algorithm = new NeedlemanWunschAffine(proteins.get(0), proteins.get(1), sm, openingGap * -1);
                ((NeedlemanWunschAffine) algorithm).setSemiGlobal(isSemiGlobal);
                ((NeedlemanWunschAffine) algorithm).setExtendingGap(-1 * extendingGap);
            } else if (chosenAlgorithm == AlgorithmType.SMITH_WATERMAN_AFFINE) {
                algorithm = new SmithWatermanAffine(proteins.get(0), proteins.get(1), sm, openingGap * -1);
                ((SmithWatermanAffine) algorithm).setExtendingGap(-1 * extendingGap);
            } else {
                algorithm = new SmithWatermanLinear(proteins.get(0), proteins.get(1), sm, openingGap * -1);
            }
            algorithm.computeAlignments();
            System.out.println("\nAlignment: \n");
            System.out.println(algorithm.getAlignment1() + " (" + proteins.get(0).getId() + ")");
            System.out.println(algorithm.getAlignment2() + " (" + proteins.get(1).getId() + ")");
            System.out.println("\nScore: " + algorithm.getScore());
            System.out.println(String.format("Identical positions: %d", algorithm.getIdenticalPositions()));
            System.out.println(String.format("Identity: %,.3f%%", algorithm.getIdentity() * 100));
            System.out.println();
        }
    }

    public static void wrongCommandExit() {
        System.out.println("\nWrong command. Usage:\n" +
                "Fasta file path: -f fastaFilePath (by default fasta.txt in directory of jar executable)\n" +
                "Alignment: -al msa (for MSA) or tsa (for two sequences alignment). By default MSA is chosen\n" +
                "Gaps: -og openingGap and -eg extendingGap (by default 10 and 1)\n" +
                "Substitution matrix: -sm smFilePath (by default blosum62.bla in directory of jar executable)\n" +
                "Algorithm for two sequences alignment: -tsa nwl (Needleman-Wunsch linear) or nwa (NW affine) or swl (Smith-Waterman linear) or swa (SM affine). By default nwl is chosen\n" +
                "Semi global for Needleman-Wunsch: -sg 0 or 1 (by default 0 which means false)\n" +
                "All pairs can be used in arbitrary sequence\n");
        System.exit(1);
    }

}
