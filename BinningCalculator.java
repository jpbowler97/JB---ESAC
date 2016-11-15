package gb.esac.myprograms;

import gb.esac.myprograms.MinMaxFreq;
import gb.esac.myprograms.KalmanFilter;
import gb.esac.myprograms.ParameterEstimator;
import gb.esac.timeseries.TimeSeries;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.timeseries.AsciiTimeSeriesFileReader;
import gb.esac.tools.LeastSquaresFitter;
import gb.esac.timeseries.TimeSeriesMaker;
import java.lang.Math;
import java.util.Arrays;


public class BinningCalculator {

	// Calculates the correct binning for a given TimeSeries
	// uses a comparison of critical frequency on processRMS values to define correct filter
	// Returns the result as {cR, nBins}
	public static double[] binningCalculator(String fileName) throws Exception {
		
	// AsciiTimeSeriesFileReader reader = new AsciiTimeSeriesFileReader();
	// // TimeSeriesFileReader reader = new TimeSeriesFileReader();

	// // extract light curve data from file as TimeSeries object - extract duration, countRate
	// TimeSeries ts = reader.readTimeSeriesFile(fileName);
	System.out.println("\nfile: " + fileName + "\n");

	TimeSeries ts = TimeSeriesMaker.makeTimeSeries(fileName);

	double leahyNorm = 2/ts.sumOfBinHeights();
	System.out.println("\nleahyNorm" + leahyNorm);

	// double[] rates = ts.getRates();
	// double ratesSum = 0;

	// for (int i = 0; i < rates.length; i++) {
	// 	ratesSum += rates[i];
	// }

	// double countRate = ratesSum/rates.length;




	double spectralIndex = ParameterEstimator.indexEstimator(ts);
	// double spectralIndex = 2.0;
	System.out.println("index = " + spectralIndex);

	double duration = ts.tStop();
	System.out.println("duration = " + duration);

	double countRate = ts.meanRate();
	// double countRate = 20;
	System.out.println("countRate = " + countRate);

	int nBins = ts.nBins;
	System.out.println("nBins (non-empty in TimeSeries) = " + nBins);

	System.out.println("\nThese are derived directly from the input file.");

	// fit to find maxFreq - make periodogram and use MinMaxFreq
	System.out.println("\nCalculating critical frequency...");

	double realMaxFreq = MinMaxFreq.maxFreq(ts);
	// double realMaxFreq = 0.5;

	// sanity check
	if (Double.isNaN(realMaxFreq)) {
		System.out.println("\nBinningCalculator - error in critical frequency calculation\n");
		System.exit(1);
	}

	// apply Kalman (look up table or calculate) to recover input index - apply Kalman filter and fit the result
	double lowBound = 0;
	double uppBound = 1.0;
	int kalIterations = 10;

	// System.out.println("critical frequency = " + realMaxFreq);	

	double processRMS = KalmanFilter.rmsIterator(ts, lowBound, uppBound, kalIterations);
	System.out.println("\nprocessRMS = " + processRMS);

	TimeSeries kalmanLC = KalmanFilter.kalmanFilter(ts, processRMS);
	
	// TEST: plots

	// kalman light curve
	System.out.println("\nWriting kalmanLC file as kalmanlc(" + fileName + ").qdp");
	kalmanLC.writeRatesAsQDP("kalmanlc(" + fileName + ").qdp");
	// System.out.println("Saving as kalmanlc.qdp");

	// original periodogram
	System.out.println("Writing unfiltered light curve as periodogram periodogram(" + fileName + ").qdp");
	PeriodogramGenerator.periodogramFileGenerator(ts, "periodogram(" + fileName + ").qdp");

	// kalman periodogram
	System.out.println("Writing kalman filtered lightcurve as periodogram kalperiodogram(" + fileName + ").qdp");
	PeriodogramGenerator.periodogramFileGenerator(kalmanLC, "kalperiodogram(" + fileName + ").qdp");

	// TEST

	double kalIndex = ParameterEstimator.indexEstimator(kalmanLC);

	// double kalIndex = 1.0;
	System.out.println("\nkalIndex = " + kalIndex); // index after filtering

	// iterate on nBins (with other input conditions fixed as calculated above) to match maxFreq - lower bound, upper bound and bisection/interpolation
	int binIterations = 10;
	int expLowBound = 5; // specifies the exponent of 2 as nBins is a power of 2
 	int expUppBound = 17;

	nBins = nBinsIterator(duration, kalIndex, countRate, realMaxFreq, binIterations, expLowBound, expUppBound);
	
	// TEST - generates periodogram based on 'correct' binning to simulate the real data set
	System.out.println("Writing as p-simulation(" + fileName + ").qdp");
	PeriodogramGenerator.periodogramFileGenerator(duration, spectralIndex, countRate, nBins, "p-simulation(" + fileName + ").qdp");
	double[] result = {countRate, (double) nBins};

	return result;
	// repeat for different countRates and plot log(nBins) against log(cR)
	}

	public static int nBinsIterator(double duration, double spectralIndex, double countRate, double realMaxFreq, int iterations, int expLowBound, int expUppBound) throws Exception {

	 	int lowBound = (int) Math.pow(2, expLowBound);
	 	int uppBound = (int) Math.pow(2, expUppBound);
	 	double maxFreq = 0;

		double[] maxFreqEsts = new double[iterations];
		double[] errors = new double[iterations];

		// Evaluate for lower and upper bounds
		System.out.println("\nIterating over nBins to match critical frequency...\n");

		// System.out.println("duration = " + duration);
		// System.out.println("index = " + spectralIndex);
		// System.out.println("countRate = " + countRate);

		System.out.println("\nnBins lowBound = " + lowBound);
		maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, lowBound);
		maxFreqEsts[0] = maxFreq;
		errors[0] = Math.abs(maxFreq - realMaxFreq);

		// record maxFreqEsts to see if wild fluctuation or convergence occurs 
		System.out.println("\nnBins uppBound = " + uppBound);
		maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, uppBound); // if uppBound is too big, causes memory error
		maxFreqEsts[1] = maxFreq;
		errors[1] = Math.abs(maxFreq - realMaxFreq);

		// this is a crude way to iterate
		for (int i = 0; i < iterations - 2; i++) {
			if (0.5*(expUppBound - expLowBound) < 0.9) { // i.e. if iterates have converged
				System.out.println("\nFinished iterating.");
				break;
			}

			if (errors[i+1] < errors[i]) {
				// expLowBound = (int) (((double) (expLowBound + expUppBound))/2.0);
				expLowBound = expLowBound + (int) (0.5*(expUppBound - expLowBound)); // slower but more accurate
				lowBound = (int) Math.pow(2, expLowBound);
				maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, lowBound);
				System.out.println("\nnBins lowBound = " + lowBound);

				maxFreqEsts[i+2] = maxFreq;
				errors[i+2] = Math.abs(maxFreq - realMaxFreq);
			}

			else {
				// expUppBound = (int) (((double) (expLowBound + expUppBound))/2.0);
				expUppBound = expUppBound - (int) (0.5*(expUppBound - expLowBound)); // slower but more accurate
				uppBound = (int) Math.pow(2, expUppBound);
				maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, uppBound);
				System.out.println("\nnBins uppBound = " + uppBound);

				maxFreqEsts[i+2] = maxFreq;
				errors[i+2] = Math.abs(maxFreq - realMaxFreq);
			}

		}

		// checks between lower and upper bound to see which is better
		maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, uppBound);
		double error1 = Math.abs(maxFreq - realMaxFreq);
		maxFreq = MinMaxFreq.maxFreq(duration, spectralIndex, countRate, lowBound);
		double error2 = Math.abs(maxFreq - realMaxFreq);

		int exponent = 0;

		if (error1 < error2) {
			exponent = expUppBound;
		}
		else { exponent = expLowBound; }

		// TEST: printing
		System.out.println("\nrealMaxFreq = " + realMaxFreq);
		System.out.println("\nmaxFreq Estimates " + Arrays.toString(maxFreqEsts) + "\n\nErrors " + Arrays.toString(errors) + "\n");
		// TEST
		System.out.println("nBins Iterator results:\nexplowBound = " + expLowBound + "\texpUppBound = " + expUppBound);

		int nBins = (int) Math.pow(2, exponent);

		return nBins;
	}

	public static void dataPlotter(double[] countRates, double[] nBins) throws Exception {
		
		// Writes plot as .qdp file
		String[] header = AsciiDataFileWriter.makeHeader("logCountRates vs logNBins", ""); // labels axes

		// Decide on fileName
		// System.out.println("Enter QDP File Name: ");
		// String fileName = reader.next();

		String fileName = "logcountRates_lognBins";


		AsciiDataFileWriter writer = new AsciiDataFileWriter(fileName + ".qdp"); // Define filename

		double[] x = new double[countRates.length];
		double[] y = new double[countRates.length];

		for (int i = 0; i < countRates.length; i++) {
			x[i] = Math.log(countRates[i]);
			y[i] = Math.log(nBins[i]);			
		}

		System.out.println("\nx values = " + Arrays.toString(x));
		System.out.println("y values = " + Arrays.toString(y));


		writer.writeData(header, x, y);

		double[] fit = LeastSquaresFitter.leastSquaresFitLine(x, y);

		// parameter sanity check
		for (int i = 0; i < fit.length; i++) {
			if (fit[i] == 0 || Double.isNaN(fit[i])) {
				System.out.println("\nError in countRates_nBins graph.\nZero or infinite valued parameter measured.");
			//	System.exit(1);
			}
		}

		System.out.println("\ngradient = " + fit[0] + "\ngradient error = " + fit[1] + "\nintercept = " + fit[2] + "\nintercept error = " + fit[3]);

		System.out.println("To look at graph, use:\nqdp '" + fileName + "'.qdp\nthen type q in command line to exit plot");

		}
}