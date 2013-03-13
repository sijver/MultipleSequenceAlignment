package core;

import core.algorithm.NeedlemanWunschAffine;
import core.algorithm.SmithWatermanAffine;
import core.io.FastaReader;
import core.io.SubstitutionMatrixReader;

import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Alignment {

    public static void main(String[] args) {
        FastaReader fr = new FastaReader("fasta.txt");
        List<Protein> proteins= fr.processFile();
//        for(Protein p : proteins){
//            System.out.println(p.toString());
//        }
        SubstitutionMatrixReader smr = new SubstitutionMatrixReader("blosum62.bla");
        SubstitutionMatrix sm = smr.processFile();
//        for(String s : sm.getAminoacids()){
//            System.out.println(s);
//        }
//        System.out.println(sm.getAminoacids().size());
//
        int[][] m = sm.getSubstitutionMatrix();
//        for(int[] d : m){
//            System.out.println(Arrays.toString(d));
//        }
        NeedlemanWunschAffine algorithm = new NeedlemanWunschAffine(proteins.get(0), proteins.get(1), sm, -10);
        algorithm.setSemiGlobal(false);
        algorithm.setExtendingGap(-1);
        algorithm.computeAlignments();
        System.out.println(algorithm.getAlignment1());
        System.out.println(algorithm.getAlignment2());
        System.out.println("Score: " + algorithm.getScore());


        m = algorithm.getSMatrix();
        for(int[] i : m){
            System.out.println(Arrays.toString(i));
        }
        System.out.println(String.format("Identical positions: %d", algorithm.getIdenticalPositions()));
        System.out.println(String.format("Identity: %,.3f%%", algorithm.getIdentity() * 100));



    }

}
