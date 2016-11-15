package gb.esac.myprograms;

import java.text.DecimalFormat;

import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;
import edu.stanford.rsl.jpop.OptimizableFunction;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.tools.LeastSquaresFitter;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.AggregatePeriodogram;
import java.lang.Math;
import java.util.Arrays;



/**
 * Fits the periodogram using the model Nx^-a + c (x is frequency) and returns the parameters N, a, c
 * Uses least squares and then minimised B-stat method
 *
 * crucial method is 'Fitter'
 * input: qdp file
 * output: estimatorParams[]
 */

// Class to fit a power law function plus constant floor
public class RedNoiseFitter implements OptimizableFunction {

    private static DecimalFormat exp = new DecimalFormat("0.000000E0");
    private static DecimalFormat num = new DecimalFormat("0.000000");
    
    double [] data;
    double [] samplingValues;
    int dimension;

    public void BFunction(double[] data, double[] samplingValues) {
	dimension = 3;
	this.data=data;
	this.samplingValues = samplingValues;
    }

    // The power law function (p = Nx^-a + c)
    // estimatorParams[0] is norm
    // estimatorParams[1] is alpha
    // estimatorParams[2] is poisson floor (constant)
    public static double powerLawFunction(double value, double[] estimatorParams) {
	return estimatorParams[0]*Math.pow(value, -estimatorParams[1]) + estimatorParams[2];
    }

    // Similarity Metric
    public double getLogLikelihoodOfModel(double[] estimatorParams) {
	double logLikelihood = 0;
	for ( int i=0; i < data.length; i++ ) {
	    logLikelihood += getLogLikelihood(powerLawFunction(samplingValues[i], estimatorParams), data[i]);
	}
	return logLikelihood;
    }

    public double getLogLikelihood(double tau, double data) {
	return -Math.log(tau) - data/tau; // tau = estimatePower
    }

    @Override
    public void setNumberOfProcessingBlocks(int number) {
	// single threaded implementation does not need this.
    }

    @Override
    public int getNumberOfProcessingBlocks() {
	return 1;
    }

    @Override
    public double evaluate(double[] x, int block) {
	double val = (-0.5 * getLogLikelihoodOfModel(x));
	if (Double.isNaN(val)) return Double.MAX_VALUE;
	return val;
    }

    // fitter methods
    public static double[] fitter(String qdpFilename) throws Exception {

		// Read restricted periodogram up to 0.01 Hz

		AsciiDataFileReader in = new AsciiDataFileReader(qdpFilename);
		
		int size = 0;
		double[] freqs = in.getDblCol(0);
		double[] powers = in.getDblCol(1);


		return fitter(freqs, powers);
	}

    public static double[] fitter(FFTPeriodogram periodogram) throws Exception {

		double[] freqs = periodogram.getFreqs();
		double[] powers = periodogram.getPowers();

		return fitter(freqs, powers);
	}

    public static double[] fitter(AveragePeriodogram periodogram) throws Exception {

		double[] freqs = periodogram.getFreqs();
		double[] powers = periodogram.getPowers();

		// System.out.println(Arrays.toString(freqs));

		return fitter(freqs, powers);
	}

	public static double[] fitter(double[] freqs, double[] powers) throws Exception {
		
		//  Read periodogram up to freqCut Hz	
		double freqCut = 0.01;
	 		
	 	//  Get initial guess from Least Squares formula
		double[] lsIndexAndNorm = fitterLS(freqs, powers, freqCut);

		double alpha = -lsIndexAndNorm[0];
		double norm = lsIndexAndNorm[1];
		if (Double.isNaN(alpha) || Double.isNaN(norm)) {
			System.out.println("\nError in RedNoiseFitter.\nLeast-squares calculation failed.");
			System.exit(1);
		}
		System.out.println("Least Squares result: \n alpha = " + alpha +"\t norm = " + norm);

		// // TEST: effectiveness of Bstat
		// double alpha = 1.8;
		// double norm = 3.4835848269932298E-5;
		// // TEST

		// double poissonFloor = 2.0;
		// double[] firstEstimators = new double[] {norm, alpha, poissonFloor};
		// String[] paramNames = new String[] {"Normalisation", "Alpha", "Poisson floor"};
	
		// int n = freqs.length;
		// double[] samplingValues = new double[n];
		// double[] data = new double[n];
		// for ( int i=0; i < n; i++ ) {
		//     samplingValues[i] = freqs[i];
		//     data[i] = powers[i];
		// }

		// TEST: Use only least-squares for speed
		double[] allEstimatorParams = new double[5]; 
		allEstimatorParams[0] = lsIndexAndNorm[1]; // norm
		allEstimatorParams[1] = -lsIndexAndNorm[0]; // alpha
		allEstimatorParams[2] = 2; // poissonFloor

		// TEST: Use only least-squares for speed

		// // Set up of function optimizer
		// FunctionOptimizer functionOptimizer = new FunctionOptimizer();
		// functionOptimizer.setDimension(3);
		// functionOptimizer.setConsoleOutput(true);
		// functionOptimizer.setOptimizationMode(OptimizationMode.Function);
		// double[] init = firstEstimators.clone();
		// functionOptimizer.setInitialX(init);

		// // Optimization
		// OptimizableFunction function = new BFunction(data, samplingValues);
		// double[] estimatorParams = functionOptimizer.optimizeFunction(function);

		// Best Model here:
		// System.out.println("Best Fit Model");
		// System.out.println("Initial B-stat: " +  function.evaluate(init, 0));
		// System.out.println("Best B-stat: " + function.evaluate(bestModel, 0));
		// System.out.println("  "+paramNames[0]+" = "+exp.format(bestModel[0]));
		// System.out.println("  "+paramNames[1]+" = "+num.format(bestModel[1]));
		// System.out.println("  "+paramNames[2]+" = "+num.format(bestModel[2]));

		// TEST: least-squares effectiveness
		// System.out.println("Least Squares result: \n alpha="+num.format(alpha)+"\t norm="+exp.format(norm));
		// TEST


		// double[] allEstimatorParams = new double[5]; // norm, alpha, poissonFloor, lsAlpha, lsNorm
		
		// allEstimatorParams[0] = estimatorParams[0]; // norm
		// allEstimatorParams[1] = estimatorParams[1]; // alpha
		// allEstimatorParams[2] = estimatorParams[2]; // poissonFloor
		allEstimatorParams[3] = alpha; // lsAlpha
		allEstimatorParams[4] = norm; // lsNorm

		// Sanity check
		System.out.println("estimatorParams = " + Arrays.toString(allEstimatorParams) + "\n(norm, index, floor, lsIndex, lsNorm)");

		for (int i = 0; i < allEstimatorParams.length; i++) {
			if (allEstimatorParams[i] <= 0) {
				System.out.println("\nRedNoiseFitter sanity check on estimated parameters failed.\nIndicates at least one parameter is negative valued.\n");
				System.exit(1);
			}
		}
		// TEST

		return allEstimatorParams;	

	      // IPlotter plotter = af.createPlotterFactory().create("Plot");
	      // plotter.createRegions(1,2,0);
	      // plotter.region(0).plot(h1d_copy);
	      // plotter.region(0).plot(fitResult1.fittedFunction());
	      // plotter.region(1).plot(h2d);
	      // plotter.show();

		//  	for (int i = 0; i < estimatorParams.length; i++){
		// 	    System.out.print(paramNames[i]+" = "+(bestModel[i])+" ");
		// 	}
		// 	System.out.println();	

	}
	public static double[] fitterLS(AveragePeriodogram periodogram) throws Exception { // least squares fits
		
		double freqCut = 0.01; // default freqCut
		return fitterLS(periodogram, freqCut);
	}

	public static double[] fitterLS(AveragePeriodogram periodogram, double freqCut) throws Exception { // least squares fits
		
		double[] freqs = periodogram.getFreqs();
		double[] powers = periodogram.getPowers();	
		return fitterLS(freqs, powers, freqCut);
	}

	public static double[] fitterLS(double[] freqs, double[] powers, double freqCut) throws Exception { // least squares fits

		//  Read periodogram up to freqCut Hz

		// System.out.println("freqs " + Arrays.toString(freqs));

	    int nValidFreqs = 0;
	    int nValidPows = 0;

	    // finds number of freqs lower than freqCut
		for (int i = 0; freqs[i] < freqCut; i++) {
			nValidFreqs++;

			// TEST
			if (powers[i] - 2 > 1E-5) { // ensures powers are above 0 so log is valid
				nValidPows++;
			}
			// TEST

			if (nValidFreqs >= freqs.length) { // prevents aioobe
				break;
			}

		}

		if (nValidFreqs < 2) { 
			System.out.println("Error in RedNoiseFitter.\nThere are no small enough frequencies to test; \nfreqCut and/or input duration are too low.");
			System.exit(1);
		}

		double[] f = new double[nValidPows];
		double[] p = new double[nValidPows];

		// System.out.println(size + " " + freqs.length);
		// System.out.println(Arrays.toString(freqs));
		// for (int i = 0; i < freqs.length; i++) {
		// 	if (freqs[i] == 0) {
		// 		System.out.println("error in RedNoiseFitter.");
		// 		System.exit(1);
		// 	}
		// }

		int i = 0;

		for (int j = 0; j < nValidFreqs; j++) {
			// f[i] = freqs[i];
			// p[i] = powers[i];

			// if (p[i] < 0) {
			// 	System.out.println("Error in RedNoiseFitter.\nNegative power read.");
			// 	System.exit(1);
			// }

			// TEST subtract expected background from powers
			if (powers[j] - 2 < 1E-5) {
				continue;
			}

			f[i] = freqs[j];
			p[i] = powers[j] - 2;
			// System.out.println(i + "\t" + p[i]);
			// TEST}

			i++;
			
			// if (size > freqs.length) { // prevents aioobe
			// 	break;
			// }
		}

		// System.out.println(Arrays.toString(f));

		if (f[nValidPows - 1] == 0) {
			System.out.println("\nError in RedNoiseFitter.\nEmpty frequency array slots.\n");
			System.exit(1);
		}

		//  Get initial guess from Least Squares formula
		double[] lsIndexAndNorm = new double[2];
		double[] result = LeastSquaresFitter.leastSquaresFitPowerLaw(f, p);

		lsIndexAndNorm[0] = result[0];
		lsIndexAndNorm[1] = result[2];

		System.out.println("Least-squares fit successful");
		
		return lsIndexAndNorm;
	}

}

