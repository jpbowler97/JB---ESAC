package gb.esac.myprograms;

import java.util.Arrays;
import java.lang.Math;

// ISSUES (doesn't finish running for some indices)
public class TestParameterEstimator {

	public static void main(String[] args) throws Exception{

		// org.apache.log4j.BasicConfigurator.configure(); // configures logger

		double duration = 10000; // duration of signal; 
		double spectralIndex = 1.5; // power law index of source1
		double countRate = 200.0; // mean event rate
		int nBins = 512; // define number of bins directly (should be pow of 2)
		double aliasingCorrectionFactor = 1;
		// double freqCut = 0.01; // used for LSfitter (default is 0.01)
		// int nFreqsPerIFS = 1; // number of frequencies per Independent Fourier Spacing, only change for oversampling
		int totalCycles = 1; // double for later calculations
		double runningAverageNorm = 0;
		double runningAverageIndex = 0;
		double runningAverageBackground = 0;
		double estimatedIndexB = 0;
		double estimatedNorm = 0;
		double estimatedBackground = 0;
		double estimatedIndexLS = 0;
		double[] averageParams = new double[5]; // norm, alpha, poissonFloor, lsAlpha, lsNorm

		// Nyquist frequency = nBins/2*duration

		// averages estimates
		for (int cycles = 0; cycles < totalCycles; cycles++) {
			double[] estimatorParams = ParameterEstimator.parameterEstimator(duration, spectralIndex, countRate, nBins, aliasingCorrectionFactor);
			System.out.println("Norm, Index, Background, lsIndex, lsNorm\n" + Arrays.toString(estimatorParams) + "\n");

			for (int i = 0; i < estimatorParams.length; i++) {
				averageParams[i] += estimatorParams[i]/totalCycles;				
			}
		}

		double expectedNorm = Math.pow((1/Math.PI), spectralIndex);

		System.out.println("countRate\t" + countRate + "\nduration\t" + duration + "\nspectralIndex\t" + spectralIndex + "\nnBins\t" + nBins + "\nexpectedNorm\t" + expectedNorm + "\n\nestimatedLSIndex\t" + averageParams[3] + "\nestimatedIndexB\t\t" + averageParams[1] + "\nestimatedNorm\t\t" + averageParams[0] + "\nestimatedBackground\t" + averageParams[2] + "\n");
	}
}