package gb.esac.myprograms;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.ChiSquare;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.binner.BinningException;
import gb.esac.eventlist.EventList;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesCombiner;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import gb.esac.tools.LeastSquaresFitter;
import gb.esac.montecarlo.RedNoiseGenerator;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.FFTPeriodogram;
import java.util.Scanner; 
import java.util.Arrays;
import org.apache.log4j.Logger;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.AggregatePeriodogram;
import gb.esac.periodogram.PeriodogramUtils;



/**
 * Generates a Red Noise periodogram and time series from event list data
 * accounts for low frequency leakage by averaging over a 128T timeSeries
 */

/**
 * Glossary:
 * Aliasing: High frequency red noiseleakage caused by not having a high enough sampling frequency
 *
 */
public class PeriodogramGenerator {

	public double acf = 1; // 128T accounts for aliasing by sampling higher frequencies; 128 was experimentally derived by G Belanger as sufficient 
	// double aliasingCorrectionFactor = 128;

	public static double getACF() {
		PeriodogramGenerator pg = new PeriodogramGenerator();
		double aliasingCorrectionFactor = pg.acf;
		return aliasingCorrectionFactor;
	}

	public static AveragePeriodogram periodogramGenerator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
		
		double aliasingCorrectionFactor = getACF();
		return periodogramGenerator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);
	}

	public static AveragePeriodogram periodogramGenerator(TimeSeries ts) throws Exception {

		AggregatePeriodogram avg = new AggregatePeriodogram();
		FFTPeriodogram p = PeriodogramMaker.makePlainFFTPeriodogram(ts, "leahy");
		avg.add(p);
		AveragePeriodogram avgPeriodogram = avg.getPeriodogram(); // .writeAsQDP("avg.qdp");
		return avgPeriodogram;
	}

	public static AveragePeriodogram periodogramGenerator(TimeSeries[] tsArray) throws Exception {

		double aliasingCorrectionFactor = getACF();		
		AggregatePeriodogram avg = new AggregatePeriodogram();


		for ( int timeSegment = 1; timeSegment <= aliasingCorrectionFactor; timeSegment++) {
			
		    FFTPeriodogram p = PeriodogramMaker.makePlainFFTPeriodogram(tsArray[timeSegment-1], "leahy");
		    // System.out.println(Arrays.toString(p.getFreqs()) + Arrays.toString(p.getPowers()));
		    avg.add(p);
		}

	    AveragePeriodogram avgPeriodogram = avg.getPeriodogram(); // .writeAsQDP("avg.qdp");
	    return avgPeriodogram;
	}

	public static AveragePeriodogram periodogramGenerator(double duration, double spectralIndex, double countRate, int nBins, int nFreqsPerIFS) throws Exception {
		double componentLossFactor = 1;
		return periodogramGenerator(duration, spectralIndex, countRate, nBins, nFreqsPerIFS, componentLossFactor);

	}

	public static AveragePeriodogram periodogramGenerator(double duration, double spectralIndex, double countRate, int nBins, int nFreqsPerIFS, double componentLossFactor) throws Exception {
		
		double aliasingCorrectionFactor = getACF();
		double effectiveDuration = aliasingCorrectionFactor*duration;

		// Remember: 
		// nuMin = 1/effectiveDuration; 
		// nuMax = freqNyq = (1/2)*samplingRate 
		// effectiveNyqFreq = 2*countRate

		TimeSeries[] tsArray = timeSeriesGenerator(duration, spectralIndex, countRate, nBins, nFreqsPerIFS, componentLossFactor);

		AveragePeriodogram avgPeriodogram = periodogramGenerator(tsArray);

		return avgPeriodogram;
	}
	

	public static AveragePeriodogram periodogramGenerator(double duration, double spectralIndex, double countRate, int nBins, double aliasingCorrectionFactor) throws Exception {

		aliasingCorrectionFactor = getACF();
		double effectiveDuration = aliasingCorrectionFactor*duration;
		int nFreqsPerIFS = 1; // number of frequencies per Independent Fourier Spacing, only change for oversampling
		AggregatePeriodogram avg = new AggregatePeriodogram();

		// Remember: 
		// nuMin = 1/effectiveDuration; 
		// nuMax = freqNyq = (1/2)*samplingRate 
		// effectiveNyqFreq = 2*countRate

		TimeSeries[] tsArray = timeSeriesGenerator(duration, spectralIndex, countRate, nBins);

		AveragePeriodogram avgPeriodogram = periodogramGenerator(tsArray);

		return avgPeriodogram;
	}

	public static String periodogramFileGenerator(AveragePeriodogram avgPeriodogram, String fileName) {
		fileName = fileName + ".qdp";
		System.out.println("\nWriting Periodogram as " + fileName + "\n");

		avgPeriodogram.writeAsQDP(fileName);

		return fileName;
	}
	public static void periodogramFileGenerator(TimeSeries ts, String fileName) throws Exception {

		AveragePeriodogram avgPeriodogram = periodogramGenerator(ts);
		// System.out.println("Saving as " + fileName);

		avgPeriodogram.writeAsQDP(fileName);	
	}

	public static String periodogramFileGenerator(double duration, double spectralIndex, double countRate, int nBins, String fileName) throws Exception {
		
		AveragePeriodogram avgPeriodogram = periodogramGenerator(duration, spectralIndex, countRate, nBins);
		System.out.println("\nSaving as " + fileName);

		avgPeriodogram.writeAsQDP(fileName);

		return fileName;		
	}

	public static String periodogramFileGenerator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
		
		AveragePeriodogram avgPeriodogram = periodogramGenerator(duration, spectralIndex, countRate, nBins);
		String fileName = Double.toString(duration) + "_" + Double.toString(spectralIndex) + "_" + Double.toString(countRate) + "_" + Integer.toString(nBins) + ".qdp";
		System.out.println("\nSaving as " + fileName);

		avgPeriodogram.writeAsQDP(fileName);

		return fileName;
	}

	// generates the periodogram as equal to the spectrum
	public static AveragePeriodogram deterministicPeriodogramGenerator(double duration, double spectralIndex, double countRate, int nBins) {

		double nuMin = 1/duration;
		double nuMax = 2*countRate;
		int nFreqsPerIFS = 1;
		// TEST
		// double nuMax = (1/2)*mean;
		// TEST
		double df = nuMin/nFreqsPerIFS;
		double nFreqs = (nuMax - nuMin)/df;

		//  Adjust nuMax to have power of 2 number of frequencies
	 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
		int nNewBins = (int) Math.pow(2, exponent);
		// if ( nFreqs != nNewBins ) {
		//     logger.warn("Number of specified frequencies ("+nFreqs+") not a power of 2. Using "+nNewBins+" instead");
		//     nFreqs = nNewBins;
		// }
		nuMax = nuMin + df*nFreqs;

		double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df);

		double[] spec = new double[frequencies.length];
		for ( int i=0; i < frequencies.length; i++ ) {
		    double omega = frequencies[i] * 2*Math.PI;
		    spec[i] = Math.pow(omega, -spectralIndex);
		}

		// generate periodogram with no variablility
		// double[] powers = spec;

		// generate periodogram with chi-squared variability
		double[] powers = new double[spec.length];
		MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
		ChiSquare chi = new ChiSquare(2.0, engine);

		for ( int i = 0; i < powers.length; i++) {
			powers[i] = (1d/2d) * spec[i] * chi.nextDouble();	
		}
		// System.out.println(Arrays.toString(powers));

		double[] errors = new double[powers.length];

		for ( int i = 0; i < errors.length; i++) {
			errors[i] = 0;	
		}

		AveragePeriodogram avgPeriodogram = new AveragePeriodogram(frequencies, powers, errors);
		return avgPeriodogram;
	}

	public static TimeSeries[] timeSeriesGenerator(double duration, double spectralIndex, double countRate, int nBins) throws Exception {
	
		int nFreqsPerIFS = 1;
		double componentLossFactor = 1;
		return timeSeriesGenerator(duration, spectralIndex, countRate, nBins, nFreqsPerIFS, componentLossFactor);
	}

	public static TimeSeries[] timeSeriesGenerator(double duration, double spectralIndex, double countRate, int nBins, int nFreqsPerIFS, double componentLossFactor) throws Exception {

		double aliasingCorrectionFactor = getACF();

		double[] times = RedNoiseGenerator.generateArrivalTimes(countRate, aliasingCorrectionFactor*duration, spectralIndex, nFreqsPerIFS, componentLossFactor);
		
		// System.out.println("\nnumber of events = " + times.length);
		// System.out.println("\nexact countRate = " + ((double)times.length)/duration);

		TimeSeries[] tsArray = new TimeSeries[(int) aliasingCorrectionFactor];
		int lastEventCounter = 0;
		int i = 0;

		for ( int timeSegment=1; timeSegment <= aliasingCorrectionFactor; timeSegment++) {
			
			int eventCounter = 0; // number of events in each TimeSeries

			// ISSUE: perhaps the try catch is now unnecessary
			try {
				while (i < times.length - 1 && times[i] <= (timeSegment*duration)) { // end point for each TimeSeries 
					eventCounter++; // nEvents in segment
					i++; // overall event number
				}
			} catch (ArrayIndexOutOfBoundsException aioobe) {
				System.out.println("i = " + i + "\n");
				System.out.println("ArrayIndexOutOfBoundsException caught in PeriodogramGenerator");
				break;
			}

			double[] actualTimes = new double[eventCounter]; // times for particular TimeSeries segment

			for(int j = 0; j < eventCounter; j++) {
				actualTimes[j] = times[j + lastEventCounter] - (timeSegment-1)*duration;
			}

			// System.out.println("\nnumber of events = " + times.length);
			// System.out.println("\nexact countRate = " + ((double)times.length)/duration);

			// ISSUES empty bins?
			EventList evlist = new EventList(actualTimes);
			tsArray[timeSegment-1] = TimeSeriesMaker.makeTimeSeries(evlist, nBins);

			// countRate = tsArray[0].meanRate();
			// System.out.println("\ncountRate = " + countRate);

			// System.out.println("nBins = " + tsArray[0].nBins);

			lastEventCounter += eventCounter;
		}

		return tsArray;
	}


}
 