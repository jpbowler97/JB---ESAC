package gb.esac.myprograms;
import java.io.File;
import java.util.Arrays;



import gb.esac.myprograms.BinningCalculator;

public class TestBinningCalculator {

	public static void main(String[] args) throws Exception {
		
		// AsciiTimeSeriesFileReader reader = new AsciiTimeSeriesFileReader;

		// TimeSeries ts = reader.readTimeSeriesFile(psd-0153950601.qdp);

		// // run through a bunch of real data files (light curves)
		// double[] countRates = new double[];
		// double[] nBins = new double[];

		// for (int i = 0; i < ; i++) {
		// 	countRates[i] = BinningCalculator.binningCalculator(ts)[0];
		// 	nBins[i] = BinningCalculator.binningCalculator(ts)[1];
		// }

		// BinningCalculator.dataPlotter(countRates, nBins);

		// ts-0153950601.qdp

		// String myDirectoryPath = "~/dev/java/gb/esac/RealDataLCs";
		// String myDirectoryPath = "./../RealDataLCs";
		String myDirectoryPath = "./../TestDataLCs";
		
		File dir = new File(myDirectoryPath);
		// System.out.println(dir.getName());
		// System.out.println(Arrays.toString(dir.list()));
		File[] directoryListing = dir.listFiles();

		if (directoryListing == null) {
			System.out.println("\nNo files could be found at " + myDirectoryPath + "\n");
		    System.exit(1);	
		}

		// System.out.println(directoryListing[0].getName());

		int i = 0;

		double[] result = new double[2];
		double[] countRates = new double[directoryListing.length]; 
		double[] nBins = new double[directoryListing.length];
		String fileName = null;

		// if (directoryListing != null) {
		    for (File lc : directoryListing) {
		        // Do something with child
		    	fileName = lc.getName();
				result = BinningCalculator.binningCalculator(fileName);// cR, nBins
				countRates[i] = result[0];
				nBins[i] = result[1];
				i++;
				System.out.println("file " + i + "\ncountRate = " + result[0] + "\tnBins = " + result[1]);
		    }
		// } else {
		//     // Handle the case where dir is not really a directory.
		//     // Checking dir.isDirectory() above would not be sufficient
		//     // to avoid race conditions with another process that deletes
		//     // directories.
		//     System.out.print("\nNo files could be found at" + myDirectoryPath + "\n");
		//     System.exit(1);
		// }

		// String[] files = ; // list of all available data sets to use 

		// String fileName = "testlc.qdp";

		// for (int i = 0; i < numberOfFiles; i ++) {
		// 	fileName = files[i];
		// 	result = BinningCalculator.binningCalculator(fileName);// cR, nBins
		// 	countRates[i] = result[0];
		// 	nBins[i] = result[1];
		// }

		BinningCalculator.dataPlotter(countRates, nBins);

	}
}