package gb.esac.myprograms;

import gb.esac.timeseries.TimeSeriesUtils;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.timeseries.TimeSeries;
import gb.esac.myprograms.PeriodogramGenerator;
import gb.esac.myprograms.LinearFitter;
import java.util.Arrays;



public class KalmanFilter {

	public static double rmsIterator(TimeSeries ts, double lowBound, double uppBound, int iterations) throws Exception {
	
		double[] indexEsts = new double[iterations];
		double[] errors = new double[iterations];

		// Evaluate for lower and upper bounds
		double[] fit = filterAndFit(ts, lowBound);

		indexEsts[0] = fit[0];
		errors[0] = fit[1];

		fit = filterAndFit(ts, uppBound);

		// record indexEsts to see if wild fluctuation or convergence occurs 
		indexEsts[1] = fit[0];
		errors[1] = fit[1];

		// this is a crude way to iterate
		for (int i = 0; i < iterations - 2; i++) {

			if (errors[i+1] < errors[i]) {
				lowBound = (lowBound + uppBound)/2.0;
				fit = filterAndFit(ts, lowBound);

				indexEsts[i+2] = fit[0];
				errors[i+2] = fit[1];
			}

			else {
				uppBound = (lowBound + uppBound)/2.0;
				fit = filterAndFit(ts, uppBound);

				indexEsts[i+2] = fit[0];
				errors[i+2] = fit[1];
			}

		}

		// TEST: printing
		// System.out.println("\nIndex Estimates " + Arrays.toString(indexEsts) + "\n\nFitting Errors " + Arrays.toString(errors));
		// TEST

		double processRMS = (uppBound + lowBound)/2.0;

		return processRMS;	

	}

	public static TimeSeries kalmanFilter(TimeSeries ts, double processRMS) throws Exception {

		TimeSeries kalmanLC = TimeSeriesUtils.kalmanFilter(ts, processRMS);		

		return kalmanLC;
	}

	// Applies the Kalman filter to an array of TimeSeries
	public static TimeSeries[] kalmanFilter(TimeSeries[] tsArray, double processRMS) throws Exception {

		TimeSeries[] kalmanLCs = new TimeSeries[tsArray.length];

		for (int i = 0; i < tsArray.length; i++) {
			TimeSeries lc = tsArray[i];
			kalmanLCs[i] = TimeSeriesUtils.kalmanFilter(lc, processRMS);		
		}

		return kalmanLCs;
	}

	public static double rmsIterator(TimeSeries[] tsArray, double lowBound, double uppBound, int iterations) throws Exception {
		
		double[] indexEsts = new double[iterations];
		double[] errors = new double[iterations];

		// Evaluate for lower and upper bounds
		double[] fit = filterAndFit(tsArray, lowBound);

		indexEsts[0] = -1*fit[0];
		errors[0] = fit[1];

		fit = filterAndFit(tsArray, uppBound);

		// record indexEsts to see if wild fluctuation or convergence occurs 
		indexEsts[1] = -1*fit[0];
		errors[1] = fit[1];

		// this is a crude way to iterate
		for (int i = 0; i < iterations - 2; i++) {

			if (errors[i+1] < errors[i]) {
				lowBound = (lowBound + uppBound)/2.0;
				fit = filterAndFit(tsArray, lowBound);

				indexEsts[i+2] = -1*fit[0];
				errors[i+2] = fit[1];
			}

			else {
				uppBound = (lowBound + uppBound)/2.0;
				fit = filterAndFit(tsArray, uppBound);

				indexEsts[i+2] = -1*fit[0];
				errors[i+2] = fit[1];
			}

		}

		// TEST: printing
		System.out.println("\nIndex Estimates " + Arrays.toString(indexEsts) + "\n\nFitting Errors " + Arrays.toString(errors));
		// TEST

		double processRMS = (uppBound + lowBound)/2.0;

		return processRMS;
	}

	public static double[] filterAndFit(TimeSeries ts, double processRMS) throws Exception {
		
		TimeSeries kalmanLC = kalmanFilter(ts, processRMS);	
		AveragePeriodogram avgPeriodogram = PeriodogramGenerator.periodogramGenerator(kalmanLC);
		double[] fit = LinearFitter.linearFitter(avgPeriodogram);

		return fit;	
	}

	public static double[] filterAndFit(TimeSeries[] tsArray, double processRMS) throws Exception {
		
		TimeSeries[] kalmanLCs = kalmanFilter(tsArray, processRMS);	
		AveragePeriodogram avgPeriodogram = PeriodogramGenerator.periodogramGenerator(kalmanLCs);
		double[] fit = LinearFitter.linearFitter(avgPeriodogram);

		return fit;
	}
}