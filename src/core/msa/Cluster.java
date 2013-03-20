package core.msa;

import core.Protein;

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

    public static List<Protein> getClusterProteins(Cluster cluster) {
        List<Protein> proteins = new LinkedList<Protein>();
        for (Object o : cluster.getClusterObjects()) {
            if (o.getClass() == Protein.class) {
                proteins.add((Protein) o);
            } else if (o.getClass() == Cluster.class) {
                proteins.addAll(Cluster.getClusterProteins((Cluster) o));
            }
        }
        return proteins;
    }

    public static String getClusterTree(Cluster cluster, int deepness){
        StringBuilder clusterTree = new StringBuilder();

        for(Object o : cluster.getClusterObjects()){
            clusterTree.append("\n");
            for(int i = 0; i < deepness; i++){
                clusterTree.append("  ");
            }
            clusterTree.append("|_");
            if(o.getClass() == Cluster.class){
                clusterTree.append(((Cluster)o).getClusterTree((Cluster)o, deepness + 1));
            } else {
                clusterTree.append(((Protein)o).getId());
            }
        }

        return clusterTree.toString();
    }
}
