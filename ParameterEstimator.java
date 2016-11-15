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
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.AggregatePeriodogram;

/**
 * Generates the estimatedSpectralIndex as  
 * a function of the periodogram and the calculated maxFreq
 *
 */

public class ParameterEstimator {

	public static double[] parameterEstimator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
		double aliasingCorrectionFactor = 128;
		return parameterEstimator(duration, spectralIndex, countRate, nBins, spectralIndex);
	}

	public static double[] parameterEstimator(double duration, double spectralIndex, double countRate, int nBins, double aliasingCorrectionFactor) throws Exception {

		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		// non-deterministic
		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);
		// deterministic
		// AveragePeriodogram periodogram = PeriodogramGenerator.deterministicPeriodogramGenerator(duration, spectralIndex, countRate, nBins);		
		// 
		return parameterEstimator(periodogram);			
	}

	public static double[] parameterEstimator(AveragePeriodogram periodogram) throws Exception {
		
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);
		System.out.println("Periodogram fitted; parameters estimated.");

		return estimatorParams;			
	}

	public static double indexEstimator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
		int aliasingCorrectionFactor = 128;
		return indexEstimator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);
	}

	public static double indexEstimator(double duration, double spectralIndex, double countRate, int nBins, double aliasingCorrectionFactor) throws Exception {

		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);

		return estimatorParams[1]; // estimated index
	}

	public static double indexEstimator(TimeSeries ts) throws Exception {

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(ts);
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);

		return estimatorParams[1]; // estimated index		
	}

	public static double normEstimator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {

		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);

		return estimatorParams[0]; // estimated index
	}

	public static double indexEstimatorB(double duration, double spectralIndex, double countRate, int nBins) throws Exception {	
		return indexEstimator(duration, spectralIndex, countRate, nBins);
	}

	public static double indexEstimatorLS(double duration, double spectralIndex, double countRate, int nBins, double freqCut) throws Exception {	
		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);
		double[] indexNorm = RedNoiseFitter.fitterLS(periodogram, freqCut);

		return indexNorm[0]; // estimated index
	}

	public static double backgroundEstimator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {

		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);
		double[] estimatorParams = RedNoiseFitter.fitter(periodogram);

		return estimatorParams[2]; // estimated index
	}

	public static double indexEstimatorLS(double duration, double spectralIndex, double countRate, int nBins) throws Exception {	

		AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins);
		double[] indexNorm = RedNoiseFitter.fitterLS(periodogram);

		return indexNorm[0]; // estimated index
	}

	// public static double[] indexEstimator(FFTPeriodogram periodogram, double spectralIndex, double countRate) throws Exception {

	// 	// org.apache.log4j.BasicConfigurator.configure(); // configures logger


	// 	double[] countEstAct = {countRate, -1.0*indexNorm[0], spectralIndex};

	// 	return countEstAct;
	// }
}
