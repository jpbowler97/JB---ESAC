package gb.esac.myprograms;

import java.text.DecimalFormat;

import edu.stanford.rsl.jpop.FunctionOptimizer;
import edu.stanford.rsl.jpop.FunctionOptimizer.OptimizationMode;
import edu.stanford.rsl.jpop.OptimizableFunction;
import gb.esac.io.AsciiDataFileReader;
import gb.esac.tools.LeastSquaresFitter;

// Class to fit a power law function plus constant floor
public class BFunction implements OptimizableFunction {

    private static DecimalFormat exp = new DecimalFormat("0.000000E0");
    private static DecimalFormat num = new DecimalFormat("0.000000");
    
    double [] data;
    double [] samplingValues;
    int dimension;

    public BFunction(double[] data, double[] samplingValues) {
	dimension = 3;
	this.data=data;
	this.samplingValues = samplingValues;
    }

    // The power law function
    // param[0] is norm
    // param[1] is alpha
    // param[2] is poisson floor (constant)
    public static double powerLawFunction(double value, double[] param) {
	return param[0]*Math.pow(value, -param[1]) + param[2];
    }

    // Similarity Metric
    public double getLogLikelihoodOfModel(double[] param) {
	double logLikelihood = 0;
	for ( int i=0; i < data.length; i++ ) {
	    logLikelihood += getLogLikelihood(powerLawFunction(samplingValues[i], param), data[i]);
	}
	return logLikelihood;
    }

    public double getLogLikelihood(double tau, double data) {
	return -Math.log(tau) - data/tau;
    }

    @Override
    public void setNumberOfProcessingBlocks(int number) {
	// single threaded implementation does not need this.
    }

    @Override
    public int getNumberOfProcessingBlocks() {
	return 1;
    }

    @Override
    public double evaluate(double[] x, int block) {
	double val = (-0.5 * getLogLikelihoodOfModel(x));
	if (Double.isNaN(val)) return Double.MAX_VALUE;
	return val;
    }

    public static void main(String args []) throws Exception {

	// Read restricted periodogram up to 0.01 Hz
    // ^^ implemented in RedNoiseFitter
	AsciiDataFileReader in = new AsciiDataFileReader("p-testperiodogram.qdp");
	double[] f = in.getDblCol(0);
	double[] p = in.getDblCol(1);

	//  Get initial guess from Least Squares formula
	double[] indexAndNorm = LeastSquaresFitter.leastSquaresFitPowerLaw(f, p);
	double alpha = -indexAndNorm[0];
	double norm = indexAndNorm[1];
	System.out.println("Least Squares result: \n alpha="+num.format(alpha)+"\t norm="+exp.format(norm));

	double poissonFloor = 2;
	double[] param = new double[] {norm, alpha, poissonFloor};
	String[] paramNames = new String[] {"Normalisation", "Alpha", "Poisson floor"};

	//  Read periodogram up to 0.1 Hz
	in = new AsciiDataFileReader("p-testperiodogram.qdp");
	f = in.getDblCol(0);
	p = in.getDblCol(1);

	int n = f.length;
	double[] samplingValues = new double[n];
	double[] data = new double[n];
	for ( int i=0; i < n; i++ ) {
	    samplingValues[i] = f[i];
	    data[i] = p[i];
	}

	// Set up of function optimizer
	FunctionOptimizer functionOptimizer = new FunctionOptimizer();
	functionOptimizer.setDimension(3);
	functionOptimizer.setConsoleOutput(true);
	functionOptimizer.setOptimizationMode(OptimizationMode.Function);
	double[] init = param.clone();
	functionOptimizer.setInitialX(init);

	// Optimization
	OptimizableFunction function = new BFunction(data, samplingValues);
	double[] bestModel = functionOptimizer.optimizeFunction(function);

	// Best Model here:
	System.out.println("Best Fit Model");
	System.out.println("Initial B-stat: " +  function.evaluate(init, 0));
	System.out.println("Best B-stat: " + function.evaluate(bestModel, 0));
	System.out.println("  "+paramNames[0]+" = "+exp.format(bestModel[0]));
	System.out.println("  "+paramNames[1]+" = "+num.format(bestModel[1]));
	System.out.println("  "+paramNames[2]+" = "+num.format(bestModel[2]));	

      // IPlotter plotter = af.createPlotterFactory().create("Plot");
      // plotter.createRegions(1,2,0);
      // plotter.region(0).plot(h1d_copy);
      // plotter.region(0).plot(fitResult1.fittedFunction());
      // plotter.region(1).plot(h2d);
      // plotter.show();

//  	for (int i = 0; i < param.length; i++){
// 	    System.out.print(paramNames[i]+" = "+(bestModel[i])+" ");
// 	}
// 	System.out.println();

    }

}
