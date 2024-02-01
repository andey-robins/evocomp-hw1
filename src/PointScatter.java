import java.io.*;

public class PointScatter extends FitnessFunction {

    String name = "Point Scatter Problem";

    public PointScatter() {
        name = "Point Scatter Problem";
    }

    public void doRawFitness(Chromo X) {
        // If we don't code enough locations for two points, the fitness
        // is not well defined, so reset
        if (Parameters.numGenes < 4) {
            System.out.println("Error: Chromosome must encode at least 2 points");
            System.exit(0);
        }

        // If we don't code an even number of locations, the fitness
        // is not well defined, so reset
        if (Parameters.numGenes % 2 != 0) {
            System.out.println("Error: Number of genes must be even");
            System.exit(0);
        }

        X.rawFitness = Double.MAX_VALUE;
        for (int i = 0; i < Parameters.numGenes; i += 2) {
            int x1 = X.getPosIntGeneValue(i);
            int y1 = X.getPosIntGeneValue(i + 1);
            for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                int x2 = X.getPosIntGeneValue(j);
                int y2 = X.getPosIntGeneValue(j + 1);
                double d = distanceBetweenPoints(x1, y1, x2, y2);
                if (d < X.rawFitness) {
                    X.rawFitness = d;
                }
            }
        }
    }

    // doPrintGenes is duplicated from the OneMax class
    public void doPrintGenes(Chromo X, FileWriter output) throws IOException {
        for (int i = 0; i < Parameters.numGenes; i++) {
            Hwrite.right(X.getGeneAlpha(i), 11, output);
        }
        output.write("   RawFitness");
        output.write("\n        ");
        for (int i = 0; i < Parameters.numGenes; i++) {
            Hwrite.right(X.getPosIntGeneValue(i), 11, output);
        }
        Hwrite.right((int) X.rawFitness, 13, output);
        output.write("\n\n");
        return;
    }

    public static double distanceBetweenPoints(int x1, int y1, int x2, int y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }
}
