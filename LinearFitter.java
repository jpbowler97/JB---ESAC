package gb.esac.myprograms;

import gb.esac.tools.LeastSquaresFitter;
import gb.esac.periodogram.AveragePeriodogram;

public class LinearFitter {
	// takes in a kalmanFiltered LC and fits in log-log space according to a line;
	// returns estimated index (i.e. slope gradient) and a measure of goodness-of-fit
	// we then iterate on the processRMS, to find the correct 
	// processRMS for best removing 'extrinsic' error
	public static double[] linearFitter(AveragePeriodogram periodogram) {
		double[] freqs = periodogram.getFreqs();
		double[] powers = periodogram.getPowers();

		double[] fit = LeastSquaresFitter.leastSquaresFitPowerLaw(freqs, powers); // fit = {index, err_index, norm, err_norm}
		double index = -1*fit[0];
		double error = fit[0];

		return fit;
	}
}