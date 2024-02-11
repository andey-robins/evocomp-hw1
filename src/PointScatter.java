import java.io.*;
import java.util.function.BiFunction;

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

        if (Parameters.dataRepresentation.equals("xycart")) {
            for (int i = 0; i < Parameters.numGenes; i += 2) {
                double x1 = mapBinary(X.getPosIntGeneValue(i), Parameters.geneSize, -1.0, 1.0);
                double y1 = mapBinary(X.getPosIntGeneValue(i + 1), Parameters.geneSize, -1.0, 1.0);

                if (isOutsideUnitCircle(x1, y1)) {
                    X.rawFitness = 0.0;
                    return;
                }

                for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                    double x2 = mapBinary(X.getPosIntGeneValue(j), Parameters.geneSize, -1.0, 1.0);
                    double y2 = mapBinary(X.getPosIntGeneValue(j + 1), Parameters.geneSize, -1.0, 1.0);

                    if (isOutsideUnitCircle(x2, y2)) {
                        X.rawFitness = 0.0;
                        return;
                    }

                    double d = distanceBetweenPoints(x1, y1, x2, y2);
                    if (d < X.rawFitness) {
                        X.rawFitness = d;
                    }
                }
            }
        } else if (Parameters.dataRepresentation.equals("thetamagpolar")) {
            for (int i = 0; i < Parameters.numGenes; i += 2) {
                double point1AngleInRadians = mapBinary(X.getPosIntGeneValue(i), Parameters.geneSize, 0.0,
                        2.0 * Math.PI);
                double point1DistanceFromOrigin = mapBinary(X.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);

                for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                    double point2AngleInRadians = mapBinary(X.getPosIntGeneValue(j), Parameters.geneSize, 0.0,
                            2.0 * Math.PI);
                    double point2DistanceFromOrigin = mapBinary(X.getPosIntGeneValue(j + 1), Parameters.geneSize, 0.0,
                            1.0);

                    double distanceBetweenPoints1And2 = polarDistance(point1AngleInRadians, point1DistanceFromOrigin,
                            point2AngleInRadians, point2DistanceFromOrigin);

                    if (distanceBetweenPoints1And2 < X.rawFitness)
                        X.rawFitness = distanceBetweenPoints1And2;
                }
            }
        } else if (Parameters.dataRepresentation.equals("degrees")) {
            for (int i = 0; i < Parameters.numGenes; i += 2) {
                // First gene represents degrees [0, 360]
                double degree = mapBinary(X.getPosIntGeneValue(i), Parameters.geneSize, 0.0, 360.0);
                // Second gene represents distance [0.0, 1.0]
                double distance = mapBinary(X.getPosIntGeneValue(i + 1), Parameters.geneSize, 0.0, 1.0);

                for (int j = i + 2; j < Parameters.numGenes; j += 2) {
                    double d = polarDistance(degree * Math.PI / 180.0, distance,
                            mapBinary(X.getPosIntGeneValue(j), Parameters.geneSize, 0.0, 360.0) * Math.PI / 180.0,
                            mapBinary(X.getPosIntGeneValue(j + 1), Parameters.geneSize, 0.0, 1.0));

                    if (d < X.rawFitness)
                        X.rawFitness = d;
                }
            }
        } else {
            System.out.println(
                    "Error: Invalid data representation parameter\nValid parameters are: `xycart`, `thetamagpolar`, and `degrees`");
            System.exit(0);
        }
    }

    public static double mapBinary(int geneValue, int bits, double min, double max) {
        // Map an integer gene value into a double of range [min, max]
        // This normalization is only safe up to ~30 bits since Java uses signed ints
        return (double) geneValue / ((1 << bits) - 1) * (max - min) + min;
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

    public static double distanceBetweenPoints(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
    }

    public static boolean isOutsideUnitCircle(double x, double y) {
        return Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2)) > 1.0;
    }

    /**
     * Perform any necessary clean up after the fitness function has finished
     * 
     * Currently only writes best fit information to best_fit.csv
     */
    public void onFinish() {
        writeChromoToCsv("best_fit.csv", "average_stats.csv", Search.bestOverAllChromo);
    }

    public static void writeChromoToCsv(String fname, String statsFilename, Chromo chromo) {
        // Define generic functions for converting from genes to points for writing
        // rather than duplicating logic for file output
        BiFunction<Integer, Integer, double[]> decodeCartesian = (g1, g2) -> {
            double x = mapBinary(g1, Parameters.geneSize, -1.0, 1.0);
            double y = mapBinary(g2, Parameters.geneSize, -1.0, 1.0);
            return new double[] { x, y };
        };
        BiFunction<Integer, Integer, double[]> decodePolar = (g1, g2) -> {
            double theta = mapBinary(g1, Parameters.geneSize, 0.0, 2.0 * Math.PI);
            double r = mapBinary(g2, Parameters.geneSize, 0.0, 1.0);
            double p1 = r * Math.cos(theta);
            double p2 = r * Math.sin(theta);
            return new double[] { p1, p2 };
        };
        BiFunction<Integer, Integer, double[]> decodeDegrees = (g1, g2) -> {
            double degree = mapBinary(g1, Parameters.geneSize, 0.0, 360.0);
            double distance = mapBinary(g2, Parameters.geneSize, 0.0, 1.0);
            double p1 = distance * Math.cos(degree * Math.PI / 180.0);
            double p2 = distance * Math.sin(degree * Math.PI / 180.0);
            return new double[] { p1, p2 };
        };

        BiFunction<Integer, Integer, double[]> decodeGenes = null;

        if (Parameters.dataRepresentation.equals("degrees")) {
            decodeGenes = decodeDegrees;
        } else if (Parameters.dataRepresentation.equals("thetamagpolar")) {
            decodeGenes = decodePolar;
        } else if (Parameters.dataRepresentation.equals("xycart")) {
            decodeGenes = decodeCartesian;
        } else {
            System.out.println(
                    "Error: Invalid data representation parameter\nValid parameters are: `xycart`, `thetamagpolar`, and `degrees`");
            System.exit(0);
        }

        try {
            FileWriter fwriter = new FileWriter(fname);

            for (int i = 0; i < Parameters.numGenes; i += 2) {
                double[] points = decodeGenes.apply(chromo.getPosIntGeneValue(i),
                        chromo.getPosIntGeneValue(i + 1));
                fwriter.write(Double.toString(points[0]));
                fwriter.write(" ");
                fwriter.write(Double.toString(points[1]));
                fwriter.write("\n");
            }

            fwriter.close();

            // Compute average average w/ std deviation
            FileWriter statsWriter = new FileWriter(statsFilename);

            if(statsWriter != null) { // Hacky way to ensure average stats aren't calculated for animations
                for(int i = 0; i < Parameters.numRuns; ++i) {
                    double averageMean = 0.0;
                    for(int j = 0; j < Parameters.generations; ++j) {
                        averageMean += Search.averageStats[i][j];
                    }

                    averageMean /= Parameters.generations;

                    double stdDeviation = 0.0;
                    for(int j = 0; j < Parameters.generations; ++j) {
                        stdDeviation += Math.pow(Search.averageStats[i][j] - averageMean, 2);
                    }

                    stdDeviation = Math.sqrt(stdDeviation / (Parameters.generations - 1));

                    
                    statsWriter.write(Double.toString(averageMean) + " " + Double.toString(stdDeviation) + "\n");
                }
            }
            
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println("Error writing best fit to file");
        }
    }
}
