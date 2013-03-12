package core.algorithm;

/**
 * Created with IntelliJ IDEA.
 */
public interface AlignmentAlgorithm {

    int getAlignment(String s1, String s2);

    void computeSMatrix();

    void computeAlignments();

    void searchMatrixMaxValue(boolean searchOverLastColumnAndRow);

    void computeIdentity();

    int getScore();

    int[][] getSMatrix();

    double getIdentity();

    int getIdenticalPositions();

    String getAlignment1();

    String getAlignment2();

}
