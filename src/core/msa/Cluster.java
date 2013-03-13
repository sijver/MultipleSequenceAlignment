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

    public static List<Integer> getAllSubclustersValues(Cluster cluster) {
        List<Integer> integers = new LinkedList<Integer>();
        for (Object o : cluster.getClusterObjects()) {
            if (o.getClass() == Integer.class) {
                integers.add((Integer) o);
            } else if (o.getClass() == Cluster.class) {
                integers.addAll(Cluster.getAllSubclustersValues((Cluster) o));
            }
        }
        return integers;
    }
}
