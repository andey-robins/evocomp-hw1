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
        } else if (Parameters.dataRepresentation.equals("thetamagpolar")){
    		for (int i = 0; i < Parameters.numGenes; i += 2) {
    			double point1AngleInRadians = mapBinary(X.getPosIntGeneValue(i), Parameters.geneSize, 0.0, 2.0 * Math.PI);
                	double point1DistanceFromOrigin = mapBinary(X.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);
        		
        		for (int j = i + 2; j < Parameters.numGenes; j += 2) {
        			double point2AngleInRadians = mapBinary(X.getPosIntGeneValue(j), Parameters.geneSize, 0.0, 2.0 * Math.PI);
                		double point2DistanceFromOrigin = mapBinary(X.getPosIntGeneValue(j + 1), Parameters.geneSize, 0.0, 1.0);
                		
            			double distanceBetweenPoints1And2 = polarDistance(point1AngleInRadians, point1DistanceFromOrigin, point2AngleInRadians, point2DistanceFromOrigin);

                    		if(distanceBetweenPoints1And2 < X.rawFitness) X.rawFitness = distanceBetweenPoints1And2;
        		}
        	}
        } else if(Parameters.dataRepresentation.equals("degrees")) {
            for(int i = 0; i < Parameters.numGenes; i += 2) {
                // First gene represents degrees [0, 360]
                double degree = mapBinary(X.getPosIntGeneValue(i), Parameters.geneSize, 0.0, 360.0);
                // Second gene represents distance [0.0, 1.0]
                double distance = mapBinary(X.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);

                for(int j = i + 2; j < Parameters.numGenes; j += 2) {
                    double d = polarDistance(degree * Math.PI / 180.0, distance,
                    mapBinary(X.getPosIntGeneValue(j), Parameters.geneSize, 0.0, 360.0) * Math.PI / 180.0,
                    mapBinary(X.getPosIntGeneValue(j + 1), Parameters.geneSize, 0.0, 1.0));

                    if(d < X.rawFitness) X.rawFitness = d;
                }
            }
        }
    	else {
    		System.out.println("Error: Invalid data representation parameter\nValid parameters are: `xycart`, `thetamagpolar`, and `degrees`");
        	System.exit(0);
    	}
}
    public static double mapBinary(int geneValue, int bits, double min, double max) {
        // Map an integer gene value into a double of range [min, max]
        // This normalization is only safe up to ~30 bits since Java uses signed ints
        return (double)geneValue / ((1 << bits) - 1) * (max - min) + min;
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

    public void onFinish() {
        // Only outputting .csv of best fit for (degrees, distance) representation currently
        if(Parameters.dataRepresentation.equals("degrees")) {
            Chromo bestChromo = Search.bestOverAllChromo;

            try {
                FileWriter bestFit = new FileWriter("best_fit.csv");

                for(int i = 0; i < Parameters.numGenes; i += 2) {
                    double degree = mapBinary(bestChromo.getPosIntGeneValue(i), Parameters.geneSize, 0.0, 360.0);
                    double distance = mapBinary(bestChromo.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);

                    bestFit.write(Double.toString(distance * Math.cos(degree * Math.PI / 180.0)));
                    bestFit.write(" ");
                    bestFit.write(Double.toString(distance * Math.sin(degree * Math.PI / 180.0)));
                    bestFit.write("\n");
                }

                bestFit.close();
            } catch(Exception e) {}
        }
        else if(Parameters.dataRepresentation.equals("thetamagpolar")) {
            Chromo bestChromo = Search.bestOverAllChromo;

            try {
                FileWriter bestFit = new FileWriter("best_fit.csv");

                for(int i = 0; i < Parameters.numGenes; i += 2) {
                    double radians = mapBinary(bestChromo.getPosIntGeneValue(i), Parameters.geneSize, 0.0, 2.0 * Math.PI);
                    double distance = mapBinary(bestChromo.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);

                    bestFit.write(Double.toString(distance * Math.cos(radians)));
                    bestFit.write(" ");
                    bestFit.write(Double.toString(distance * Math.sin(radians)));
                    bestFit.write("\n");
                }

                bestFit.close();
            } catch(Exception e) {}
        }
    }
}
