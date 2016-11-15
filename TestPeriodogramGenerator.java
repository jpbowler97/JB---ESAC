package gb.esac.myprograms;

import java.util.Scanner; 
import gb.esac.periodogram.FFTPeriodogram;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.periodogram.AggregatePeriodogram;
import gb.esac.timeseries.TimeSeries;


public class TestPeriodogramGenerator {

	public static void main(String[] args) throws Exception {

		// Scanner reader = new Scanner(System.in); // records user input
		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		double countRate = 20.0; // mean event rate
		double duration = 10000; // duration of signal; 
		double spectralIndex = 2.0; // power law index of source1
		int nBins = 1024; // define number of bins directly (should be pow of 2)
		double aliasingCorrectionFactor = 1; // ISSUE: setting to 1 causes ioobe 
		int deterministic = 0; // flag for deterministic periodogram
		
		AveragePeriodogram avgPeriodogram = new AveragePeriodogram();

		// Nyquist frequency = nBins/2*duration
		// double spectralIndex2 = 2; // power law index of source2
		// double nuBreak1 = 0.1; // frequency threshold between source1 and source2
		// int nFreqsPerIFS = 1; // number of frequencies per Independent Fourier Spacing, only change for oversampling
		
		// Could add to implement below (easily done by seeing deterministicPeriodogram method)
		// Produces a deterministic periodogram i.e. the spectrum
		// if (deterministic == 2) {
		// 	avgPeriodogram = PeriodogramGenerator.superDeterministicPeriodogramGenerator(duration, spectralIndex, countRate, nBins);
		// }

		// Produces a deterministic periodogram i.e. (1/2)*spectrum*chi-2
		if (deterministic == 1) {
			avgPeriodogram = PeriodogramGenerator.deterministicPeriodogramGenerator(duration, spectralIndex, countRate, nBins);
		}
		else {			// Simulates real periodogram
			avgPeriodogram = PeriodogramGenerator.periodogramGenerator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);	
		}

		// name file yourself
		// System.out.println("Enter QDP File Name: ");
		// String fileName = reader.next();
		String fileName = "testperiodogram";

		// // automatically names file
		// String fileName = "_d" + Double.toString(duration) + "_si" + Double.toString(spectralIndex) + "_cr" + Double.toString(countRate) + "_b" + Integer.toString(nBins);
		
		if (deterministic == 1) {
			fileName += "deterministic";
		}

		/**
		 * Writes periodogram as a .qdp file
		 */

		String pFileName = "p-" + fileName;
		avgPeriodogram.writeAsQDP(pFileName + ".qdp");

		System.out.println("\nThe periodogram file name is " + pFileName + "\nTo look at file, use:\ncd ~/dev/java/gb/esac/myprograms; qdp p-" + fileName + ".qdp\n"); // \nUse the /GIF PGPLOT type, then type EXIT and launch through finder");
		
		// to look at file:
		// qdp p-fileName.qdp
		// 

		// if qdp command isn't recongised:
		// export HEADAS=/Heasoft/heasoft-6.19/x86_64-apple-darwin15.4.0
  		// . $HEADAS/headas-init.sh

	// // TEST TimeSeries Generator
	// TimeSeries[] tsArray = PeriodogramGenerator.timeSeriesGenerator(10000, 2.0, 150, 256); // avoid empty bins
	// String fileName = "test";
	// tsArray[0].writeRatesAsQDP(fileName + "lc21.qdp");

	// // TEST PeriodogramGenerator
	// PeriodogramGenerator.periodogramFileGenerator(tsArray[0], fileName + "periodogram.qdp");

	}
}
