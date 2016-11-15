package gb.esac.myprograms;

import cern.colt.list.DoubleArrayList;
import gb.esac.binner.BinningException;
import gb.esac.eventlist.EventList;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import gb.esac.tools.LeastSquaresFitter;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.FFTPeriodogram;
import java.util.Scanner; 
import java.util.Arrays;
import java.lang.Math;
import org.apache.log4j.Logger;
import gb.esac.myprograms.RedNoiseFitter;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.AggregatePeriodogram;



public class MinMaxFreq {

	// private static Logger logger  = Logger.getLogger(MinMaxFreq.class);

	// public static double maxFreq(double countRate, double duration, double spectralIndex, int nBins) throws Exception {

	// 	// org.apache.log4j.BasicConfigurator.configure(); // configures logger

	// 	FFTPeriodogram fft = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);
	// 	double[] powArray = fft.getPowers();
	// 	double[] freqArray = fft.getFreqs();

	// 	int count = 0;

	// 	while(powArray[count] > 2) {
	// 		count++;
	// 	}

	// 	return freqArray[count];
	// }

	// public static double minFreq(double countRate, double duration, double spectralIndex, int nBins) throws Exception {
		
	// 	FFTPeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);

	// 	return minFreq(periodogram, spectralIndex);

	// }
	// public static double minFreq(FFTPeriodogram periodogram, double spectralIndex) throws Exception {
		
	// 	// org.apache.log4j.BasicConfigurator.configure(); // configures logger

	// 	double[] powers = periodogram.getPowers();
	// 	double[] freqs = periodogram.getFreqs();

	// 	int count = 0;

	// 	try {	
	// 		while(powers[count] < Math.pow(1/(2*Math.PI*freqs[count]), spectralIndex)) { // compares power values with perfect power law line
	// 			count++;
	// 		}
	// 	} catch (IndexOutOfBoundsException ioob) {
	// 		System.out.println("All powers are below power-law line.\nNo minimum frequency.");
	// 		System.out.println("Program exit");
	// 		throw ioob;
	// 	}	
		
	// 	System.out.println("minFreq = " + freqs[count]);
	// 	return freqs[count]; // first frequency which has power above the power law line
	// }

	public static double maxFreq(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
		double nFreqsPerIFS = 1;
		return maxFreq(duration, spectralIndex, countRate, nBins, nFreqsPerIFS);	
	}

	public static double maxFreq(double duration, double spectralIndex, double countRate, int nBins, double componentLossFactor) throws Exception {
		int nFreqsPerIFS = 1;
		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins, nFreqsPerIFS, componentLossFactor);

		return maxFreq(periodogram);
	}

	public static double maxFreq(TimeSeries ts) throws Exception {

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(ts);
		return maxFreq(periodogram);
	}

	public static double maxFreq(FFTPeriodogram periodogram) throws Exception {

		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);
		double maxFreq = Math.pow((estimatorParams[0]/estimatorParams[2]), 1/estimatorParams[1]); // x = (N/c)^(1/a)

		return maxFreq; //spectral index

	}
	
	public static double maxFreq(AveragePeriodogram periodogram, String fileLabel) throws Exception {
		
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);
		double maxFreq = Math.pow((estimatorParams[0]/estimatorParams[2]), 1/estimatorParams[1]); // x = (N/c)^(1/a)
		System.out.println("\ncritical frequency = " + maxFreq);

		// TEST create periodogram file
		PeriodogramGenerator.periodogramFileGenerator(periodogram, "p-" + fileLabel + "grpF_" + maxFreq);
		// TEST

		return maxFreq; //spectral index		
	}

	public static double maxFreq(AveragePeriodogram periodogram) throws Exception {

		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);
		double maxFreq = Math.pow((estimatorParams[0]/estimatorParams[2]), 1/estimatorParams[1]); // x = (N/c)^(1/a)
		System.out.println("\ncritical frequency = " + maxFreq);

		// TEST create periodogram file
		PeriodogramGenerator.periodogramFileGenerator(periodogram, "p-grpF_" + maxFreq);
		// TEST

		return maxFreq; //spectral index

	}
}