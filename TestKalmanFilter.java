package gb.esac.myprograms;

import gb.esac.myprograms.KalmanFilter;
import gb.esac.myprograms.PeriodogramGenerator;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.timeseries.TimeSeries;

public class TestKalmanFilter{
	// Constructs a Kalman filtered LC based on specified input conditions;
	public static void main(String[] args) throws Exception {
		double countRate = 2.0; 
		double duration = 2000; 
		double spectralIndex = 4.0; 
		int nBins = 1024; 
		double lowBound = 0.1;
		double uppBound = 2.0;
		double processRMS = 0.6;
		int iterations = 15;

		TimeSeries[] tsArray = PeriodogramGenerator.timeSeriesGenerator(duration, spectralIndex, countRate, nBins);

		// // Applies filter and generates Periodogram, which we can then compare to the original
		// AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(tsArray);
		// String originalPeriodogramFile = PeriodogramGenerator.periodogramFileGenerator(periodogram, "originalPeriodogramFile");

		// TimeSeries[] kalmanLCs = KalmanFilter.kalmanFilter(tsArray, processRMS);
		// AveragePeriodogram kalPeriodogram = PeriodogramGenerator.periodogramGenerator(kalmanLCs);
		// String kalmanPeriodogramFile = PeriodogramGenerator.periodogramFileGenerator(kalPeriodogram, "kalmanPeriodogramFile");

		// System.out.println("\nThe original periodogram file name is " + originalPeriodogramFile + "\nTo look at file, use:\ncd ~/dev/java/gb/esac/myprograms; qdp " + originalPeriodogramFile); // \nUse the /GIF PGPLOT type, then type EXIT and launch through finder");
		// System.out.println("\n\nThe kalman periodogram file name is " + kalmanPeriodogramFile + "\nTo look at file, use:\ncd ~/dev/java/gb/esac/myprograms; qdp " + kalmanPeriodogramFile + "\n"); // \nUse the /GIF PGPLOT type, then type EXIT and launch through finder");

		// // Iterates to find optimal processRMS to remove white noise
		processRMS = KalmanFilter.rmsIterator(tsArray, lowBound, uppBound, iterations);
		System.out.println("\nOptimal processRMS = " + processRMS);
	}
}