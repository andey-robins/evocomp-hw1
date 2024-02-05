import java.io.*;

public class PointScatter extends FitnessFunction {

    String name = "Point Scatter Problem";

    public PointScatter() {
        name = "Point Scatter Problem";
    }

    public void doRawFitness(Chromo X) {
        if (Parameters.numGenes < 4) {
            System.out.println("Error: Chromosome must encode at least 2 points");
            System.exit(0);
        }

        if (Parameters.numGenes % 2 != 0) {
            System.out.println("Error: Number of genes must be even");
            System.exit(0);
        }

        X.rawFitness = Double.MAX_VALUE;
        
        
        if (Parameters.dataRepresentation.equals("xycart")){
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
        else if (Parameters.dataRepresentation.equals("thetamagpolar")){
    		for (int i = 0; i < Parameters.numGenes; i += 2) {
        		double theta1 = X.getPosIntGeneValue(i) * (2 * Math.PI / 360);
        		double r1 = X.getPosIntGeneValue(i + 1);
        		for (int j = i + 2; j < Parameters.numGenes; j += 2) {
            			double theta2 = X.getPosIntGeneValue(j) * (2 * Math.PI / 360);
            			double r2 = X.getPosIntGeneValue(j + 1);

            			double d = polarDistance(theta1, r1, theta2, r2);
            			if (d < X.rawFitness) {
                			X.rawFitness = d;
            			}
        		}
        	}
    	}
    	else {
    		System.out.println("Error: Invalid data representation parameter\nValid parameters are: `xycart` and `thetamagpolar`");
        	System.exit(0);
    	}
}

public static double polarDistance(double theta1, double r1, double theta2, double r2) {
    double x1 = r1 * Math.cos(theta1);
    double y1 = r1 * Math.sin(theta1);
    double x2 = r2 * Math.cos(theta2);
    double y2 = r2 * Math.sin(theta2);

    return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
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
