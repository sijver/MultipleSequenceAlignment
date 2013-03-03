package core;

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
        SubstitutionMatrixReader smr = new SubstitutionMatrixReader("pam250.bla");
        SubstitutionMatrix sm = smr.processFile();
//        for(String s : sm.getAminoacids()){
//            System.out.println(s);
//        }
//        System.out.println(sm.getAminoacids().size());
//
//        int[][] m = sm.getSubstitutionMatrix();
//        for(int[] d : m){
//            System.out.println(Arrays.toString(d));
//        }
        NeedlemanWunsch nw = new NeedlemanWunsch(proteins.get(0), proteins.get(1), sm, -14);
        nw.computeAlignments();
        System.out.println(nw.getAlignment1());
        System.out.println(nw.getAlignment2());
        System.out.println("Score: " + nw.getScore());

        int[][] m = nw.getfMatrix();
        for(int[] i : m){
            System.out.println(Arrays.toString(i));
        }
        System.out.println(String.format("Identical positions: %d", nw.getIdenticalPositions()));
        System.out.println(String.format("Identity: %,.3f%%", nw.getIdentity() * 100));



    }

}
