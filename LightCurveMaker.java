package gb.esac.myprograms;

import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister64;
import gb.esac.binner.BinningException;
import gb.esac.binner.BinningUtils;
import gb.esac.io.AsciiDataFileWriter;
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.PeriodogramMaker;
import gb.esac.periodogram.PeriodogramUtils;
import gb.esac.periodogram.WindowFunctionException;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.tools.BasicStats;
import gb.esac.tools.Complex;
import gb.esac.tools.FFT;
import java.text.DecimalFormat;
import java.util.Date;
import org.apache.log4j.Logger;
import java.util.Arrays;
import java.io.IOException;
import gb.esac.montecarlo.TimmerKonig;


public class LightCurveMaker {

    static Logger logger = Logger.getLogger(LightCurveMaker.class);

	public static void main(String[] args) throws Exception {
	  
		double mean = 2.0; // mean event rate
		double duration = 6000; // duration of signal; 
		double alpha = 64.0; // power law index of source1
		int nBins = 1024; // define number of bins directly (should be pow of 2)
		double aliasingCorrectionFactor = 128;
		double nFreqsPerIFS = 1;

		double nuMin = 1/duration;
		double nuMax = 2*mean;
		// TEST
		// double nuMax = (1/2)*mean;
		// TEST
		double df = nuMin/nFreqsPerIFS;
		double nFreqs = (nuMax - nuMin)/df;

		// Adjust nuMax to have power of 2 number of frequencies
	 	double exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
		int nNewBins = (int) Math.pow(2, exponent);
		if ( nFreqs != nNewBins ) {
		    logger.warn("Number of specified frequencies ("+nFreqs+") not a power of 2. Using "+nNewBins+" instead");
		    nFreqs = nNewBins;
		}
		nuMax = nuMin + df*nFreqs;
		
		// Generate Fourier components, get the rates, and scale them to the specified mean rate
		// Get components non-deterministically
		// Complex[] fourierComp = getFourierComponentsForFrequencies(alpha, nuMin, nuMax, df); 
		// Get components deterministically - takes much longer
		// Complex[] fourierComp = generateComponentsDeterministically(alpha, nuMin, nuMax, df); 		
		
		// TEST: take differences of fourier components between deterministic and non-deterministic
		Complex[] fourierComp1 = TimmerKonig.generateComponentsDeterministically(alpha, nuMin, nuMax, df); 
		double df1 = df;





		// nuMin = 1/(aliasingCorrectionFactor*duration);
		// nuMax = 2*mean;
		// // TEST
		// // double nuMax = (1/2)*mean;
		// // TEST
		// df = nuMin/nFreqsPerIFS;
		// nFreqs = (nuMax - nuMin)/df;

		// //  Adjust nuMax to have power of 2 number of frequencies
	 // 	exponent = Math.floor(Math.log10(nFreqs)/Math.log10(2));
		// nNewBins = (int) Math.pow(2, exponent);
		// if ( nFreqs != nNewBins ) {
		//     logger.warn("Number of specified frequencies ("+nFreqs+") not a power of 2. Using "+nNewBins+" instead");
		//     nFreqs = nNewBins;
		// }
		// nuMax = nuMin + df*nFreqs;

		Complex[] fourierComp2 = TimmerKonig.getFourierComponentsForFrequencies(alpha, nuMin, nuMax, df); 		
		double df2 = df;




		// Complex[] fourierComp = new Complex[fourierComp2.length];
		Complex[] fourierComp = new Complex[fourierComp1.length];

		// System.out.println(fourierComp1.length - fourierComp2.length);
		// System.out.println(fourierComp1[10]);

		for ( int i=0; i < fourierComp.length; i++ ) {
			// fourierComp[i] = fourierComp2[i].minus(fourierComp1[i/128]); // nonDet - det
			fourierComp[i] = fourierComp1[i];
			// fourierComp[i] = fourierComp2[i];			
		}
		// TEST

		double[] rates = TimmerKonig.getRatesFromFourierComponents(fourierComp);
		double scalingFactor = mean/BasicStats.getMean(rates);
		for ( int i=0; i < rates.length; i++ ) {
		    rates[i] *= scalingFactor;
		}
		// System.out.println(Arrays.toString(rates));

		// TEST: Plots Light Curve
		double[] x_Axis = new double[rates.length];
		double[] y_Axis = new double[rates.length];

		for ( int i=0; i < rates.length; i++ ) {
		    // x_Axis[i] = (aliasingCorrectionFactor*duration/rates.length)*i; // non-det time
		    x_Axis[i] = (duration/rates.length)*i; // det time 
		    y_Axis[i] = rates[i]; // rates
		}	

		// System.out.println("Hello");
		String fileName = "hist" + "_d" + Double.toString(duration) + "_si" + Double.toString(alpha) + "_cr" + Double.toString(mean);
		String[] header = AsciiDataFileWriter.makeHeader("Det vs Non-Det discrepancy", ""); // labels axes
		AsciiDataFileWriter writer = new AsciiDataFileWriter(fileName + ".qdp"); // Define filename
		writer.writeData(header, x_Axis, y_Axis);
		// TEST

		System.out.println("To view light curve, use:\nqdp '" + fileName + "'.qdp\nthen type q in command line to exit plot");

	}
}