package core.msa;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Cluster {

    private List<Object> clusterObjects;

    public Cluster() {
        clusterObjects = new LinkedList<Object>();
    }

    public void addObjectToCluster(Object object) {
        clusterObjects.add(object);
    }

    public List<Object> getClusterObjects() {
        return clusterObjects;
    }

    public static List<Integer> getClusterProteins(Cluster cluster) {
        List<Integer> integers = new LinkedList<Integer>();
        for (Object o : cluster.getClusterObjects()) {
            if (o.getClass() == Integer.class) {
                integers.add((Integer) o);
            } else if (o.getClass() == Cluster.class) {
                integers.addAll(Cluster.getClusterProteins((Cluster) o));
            }
        }
        return integers;
    }

    public static String getClusterTree(Cluster cluster, int deepness){
        StringBuilder clusterTree = new StringBuilder();

        clusterTree.append("\n");
        for(int i = 0; i < deepness; i++){
            clusterTree.append("  ");
        }
        clusterTree.append("|_");
        for(Object o : cluster.getClusterObjects()){
            if(o.getClass() == Cluster.class){
                clusterTree.append(((Cluster)o).getClusterTree((Cluster)o, deepness + 1));
            } else {
                clusterTree.append(o);
            }
        }

        return clusterTree.toString();
    }
}
