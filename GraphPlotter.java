package gb.esac.myprograms;

import gb.esac.io.AsciiDataFileWriter;
import java.util.Scanner; 
import java.util.Arrays;
import gb.esac.periodogram.AveragePeriodogram;
import gb.esac.tools.LeastSquaresFitter;
import java.io.PrintWriter;





/**
 * Plots a graph for various input values
 * and the resulting maximum frequency (according to MinMaxFreq.maxFreq() method)
 *
 * ISSUE: multiply powers by 2 to simulate Leahy normalisation?
 */

public class GraphPlotter {

	public static void main(String[] args) throws Exception {

		Scanner reader = new Scanner(System.in); // records user input

		double[] duration = {200, 500, 1000, 2000, 3000, 4000, 10000}; // less than or equal to 100 can cause error with least squares as there are no sufficiently low frequencies (see RedNoiseFitter)
		double[] indices = {1.0, 1.3, 1.5, 1.8, 2.1, 2.5};
		double[] countRate = {2.0, 4.0, 8.0, 16.0, 20.0, 32.0, 64.0, 128.0, 256.0};
		int[] nBins = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 262144};
		double[] correctionFactors = {1, 2, 4, 32, 64, 128, 256, 512, 1024}; // derive periodogram from Timeseries of length correctionFactors*T (where T is the duration you wish to simulate) - 128 eliminates aliasing?
		int[] nFreqsPerIFS = {1, 2, 4, 8};
		double[] componentLossFactor = {0.1, 10.0, 100.0, 1000.0, 100000.0};

		// double[] correctionFactors = new int[duration.length];
		// for (int i = 0; i < duration.length; i++) {
		// 	correctionFactors[i] = 250000/duration[i]; // generates all signals based on a constant duration of 250,000
		// }

		double[] x = null;
		double[] y = null;
		int[] xInt = null;
		double[] maxFreqs = null;
		double[] estimatedIndices = null;
		double[] estimatedIndexErrors = null;
		double[] maxPower = null;
		double runningAverageIndex = 0;
		double runningAveragePower = 0;
		double runningAverageFreq = 0;
		double nuMin = 0;
		double nuMax = 0;
		double df = 0;

		// Nyquist frequency = nBins/2*duration

		// choose default inputs
		double totalCycles = 1; // number of simulated data sets for which each y coordinate is averaged (double for division later)
		double d = duration[6];
		double si = indices[3];
		double cr = countRate[4];
		int b = nBins[9];
		double cf = correctionFactors[0];
		double deterministic = 0; // set 1 for deterministic 0 otherwise
		int nf = nFreqsPerIFS[0];
		double clf = componentLossFactor[0];

		// enter plot type here
		String plot = "nFreqsPerIFS-maxFreq"; // define independent variable (see cases below)

		switch (plot) { // computes y coordinates based on ^^; if the choice was invalid, it will say so.
			case "index-maxFreq": 
				x = indices;
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, indices[i], cr, b);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;	
					System.out.println("index: " + indices[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
				}
				y = maxFreqs;
				break;

			case "componentLossFactor-maxFreq": 
				x = componentLossFactor;
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, cr, b, componentLossFactor[i]);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;	
					System.out.println("componentLossFactor: " + componentLossFactor[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
				}
				y = maxFreqs;
				break;

			case "log10grpFreq-log10componentLossFactor": 
				x = componentLossFactor;
				double[] y2 = new double[x.length];
				double[] x2 = new double[x.length];
				double[] grpFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, cr, b, componentLossFactor[i]);
						// will also generate Time Domain Signals and periodograms
					}
					grpFreqs[i] = runningAverageIndex/totalCycles;	
					System.out.println("componentLossFactor: " + componentLossFactor[i] + "\t\tgrpFreq: " + grpFreqs[i]);
					runningAverageIndex = 0;
					y2[i] = Math.log10( grpFreqs[i] );
					x2[i] = Math.log10( x[i] );
				}
				x = x2;
				y = y2;

				System.out.println("\nLeast-squares fitting according to a power law...");
				double[] fit = LeastSquaresFitter.leastSquaresFitLine(x2, y2);	
				System.out.println("\nGradient = " + fit[0] + "\nGradient error = " + fit[1] + "\nIntercept = " + fit[2] + "\nIntercept error = " + fit[3]);
	

				break;

			case "log10grpFreq-log10nFourCpts": 
				x = componentLossFactor;
				y2 = new double[x.length];
				x2 = new double[x.length];
				grpFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageFreq += MinMaxFreq.maxFreq(d, si, cr, b, componentLossFactor[i]);
						// will also generate Time Domain Signals and periodograms
					}
					grpFreqs[i] = runningAverageFreq/totalCycles;	
					System.out.println("componentLossFactor: " + componentLossFactor[i] + "\t\tgrpFreq: " + grpFreqs[i]);
					runningAverageFreq = 0;

					nuMin = 1.0/d;
					nuMax = cr;
					df = nuMin * componentLossFactor[i];
					x[i] = (nuMax - nuMin)/df; // nFreqs

					y2[i] = Math.log10( grpFreqs[i] );
					x2[i] = Math.log10( x[i] );
				}
				x = x2;
				y = y2;

				System.out.println("\nLeast-squares fitting according to a power law...");
				fit = LeastSquaresFitter.leastSquaresFitLine(x2, y2);	
				System.out.println("\nGradient = " + fit[0] + "\nGradient error = " + fit[1] + "\nIntercept = " + fit[2] + "\nIntercept error = " + fit[3]);
	

				break;

			case "3D-x_index_y_log10clf_z_log10grpFreq": 
				
				x = indices;
				// double[] indexEsts = new double[x.length];
				y = componentLossFactor;
				double[][] z2 = new double[x.length][y.length];
				y2 = new double[y.length]; // plots clf
				// double[] y3 = new double[y.length]; // plots nFourierCpts 
				double[][] groupFreqs = new double[x.length][y.length];

				for (int i = 0; i < x.length; i++) {
					for (int j = 0; j < y.length; j++) {
						for (int cycles = 0; cycles < totalCycles; cycles++) {
							AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(d, indices[i], cr, b, nf, componentLossFactor[j]);
							String fileLabel = "_si" + indices[i] + "_clf" + componentLossFactor[j] + "_nf" + nf;
							runningAverageFreq += MinMaxFreq.maxFreq(periodogram, fileLabel);
							// will also generate Time Domain Signals and periodograms

							// alternatively plot index as the average predicted index
							// runningAverageIndex += RedNoiseFitter.fitter(periodogram)[0];
						}
						// indexEsts[i] = runningAverageIndex/totalCycles;
						groupFreqs[i][j] = runningAverageFreq/totalCycles;	
						System.out.println("index: "  + indices[i] + "\tcomponentLossFactor: " + componentLossFactor[j] + "\t\tgrpFreq: " + groupFreqs[i][j]);
						// alternatively plot index as the average predicted index
						runningAverageFreq = 0;
						runningAverageIndex = 0;
						z2[i][j] = Math.log10( groupFreqs[i][j] );
						y2[j] = Math.log10( y[j] );

						// double df = (nuMin/nFreqsPerIFS) * componentLossFactor;
						// double nFreqs = (nuMax - nuMin)/df;	
						// nuMin = 1/duration;
						// nuMax = cr;
						// df = nuMin * componentLossFactor[j];
						// y3[j] = (nuMax - nuMin)/df; // nFreqs
					}
				}
				// x = indexEsts;
				y = y2;
				// y = y3;
				double[][] z = z2;

				System.out.println("\n" + Arrays.toString(x) + "\n\n" + Arrays.toString(y) + "\n\n" + Arrays.toString(z[0]) + "\n\n" + Arrays.toString(z[1]) + "\n\n" + Arrays.toString(z[2]) + "\n\n" + Arrays.toString(z[3]));
				
				PrintWriter writer = new PrintWriter("the-file-name.txt", "UTF-8");
				// writer.println(x[0] + "\t" + y[0] + "\t" + z[0][0]);
				// writer.println(x[0] + "\t" + y[1] + "\t" + z[0][1]);
				for (int i = 0; i < x.length; i++) {
					for (int j = 0; j < y.length; j++) {
						writer.println(x[i] + "\t" + y[j] + "\t" + z[i][j] + "\t" + 0.0);
					}
					writer.println("");
				}

				writer.close();

				System.exit(0);

				break;

			case "3D-x_index_y_log10nFourCpts_z_log10grpFreq": 
				
				x = indices;
				// double[] indexEsts = new double[x.length];
				y = componentLossFactor;
				z2 = new double[x.length][y.length];
				y2 = new double[y.length]; // plots clf
				double[] y3 = new double[y.length]; // plots nFourierCpts 
				groupFreqs = new double[x.length][y.length];

				for (int i = 0; i < x.length; i++) {
					for (int j = 0; j < y.length; j++) {
						for (int cycles = 0; cycles < totalCycles; cycles++) {
							AveragePeriodogram periodogram = PeriodogramGenerator.periodogramGenerator(d, indices[i], cr, b, nf, componentLossFactor[j]);
							String fileLabel = "_si" + indices[i] + "_clf" + componentLossFactor[j];
							runningAverageFreq += MinMaxFreq.maxFreq(periodogram, fileLabel);
							// will also generate Time Domain Signals and periodograms

							// alternatively plot index as the average predicted index
							// runningAverageIndex += RedNoiseFitter.fitter(periodogram)[0];
						}
						// indexEsts[i] = runningAverageIndex/totalCycles;
						groupFreqs[i][j] = runningAverageFreq/totalCycles;	
						System.out.println("index: "  + indices[i] + "\tcomponentLossFactor: " + componentLossFactor[j] + "\t\tgrpFreq: " + groupFreqs[i][j]);
						// alternatively plot index as the average predicted index
						runningAverageFreq = 0;
						runningAverageIndex = 0;
						z2[i][j] = Math.log10( groupFreqs[i][j] );

						// double df = (nuMin/nFreqsPerIFS) * componentLossFactor;
						// double nFreqs = (nuMax - nuMin)/df;	
						nuMin = 1.0/d;
						nuMax = cr;
						df = nuMin * componentLossFactor[j];
						y3[j] = (nuMax - nuMin)/df; // nFreqs
						y2[j] = Math.log10( y3[j] );
					}
				}
				// x = indexEsts;
				y = y2;
				y = y3;
				z = z2;

				System.out.println("\n" + Arrays.toString(x) + "\n\n" + Arrays.toString(y) + "\n\n" + Arrays.toString(z[0]) + "\n\n" + Arrays.toString(z[1]) + "\n\n" + Arrays.toString(z[2]) + "\n\n" + Arrays.toString(z[3]));
				
				writer = new PrintWriter("the-file-name.txt", "UTF-8");
				// writer.println(x[0] + "\t" + y[0] + "\t" + z[0][0]);
				// writer.println(x[0] + "\t" + y[1] + "\t" + z[0][1]);
				for (int i = 0; i < x.length; i++) {
					for (int j = 0; j < y.length; j++) {
						writer.println(x[i] + "\t" + y[j] + "\t" + z[i][j] + "\t" + 0.0);
					}
					writer.println("");
				}

				writer.close();

				System.exit(0);

				break;

			case "log10componentLossFactor-log10FreqRatio": 
				x = componentLossFactor;
				y2 = new double[x.length];
				x2 = new double[x.length];
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, cr, b, componentLossFactor[i]);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;	
					System.out.println("componentLossFactor: " + componentLossFactor[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
					y2[i] = Math.log10( maxFreqs[i]/cr );
					x2[i] = Math.log10( x[i] );
				}
				x = x2;
				y = y2;

				System.out.println("\nLeast-squares fitting according to a power law...");
				fit = LeastSquaresFitter.leastSquaresFitLine(x2, y2);	
				System.out.println("\nGradient = " + fit[0] + "\nGradient error = " + fit[1] + "\nIntercept = " + fit[2] + "\nIntercept error = " + fit[3]);
	

				break;

			case "nBins-maxFreq": 
				xInt = nBins; 
				x = new double[xInt.length];
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, cr, nBins[i]);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;	
					x[i] = (double) xInt[i]; // converts int to double to agree with other cases
					System.out.println("nBins: " + nBins[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
				}	
				y = maxFreqs;		
				break;

			case "nFreqsPerIFS-maxFreq": 
				xInt = nFreqsPerIFS;
				x = new double[xInt.length];
				maxFreqs = new double[x.length];
				String fileLabel = null;
				AveragePeriodogram periodogram = null;

				// In Progress
				double nuMin = 1/(duration); // should be 128*duration
				double nuMax = mean;
				double df = (nuMin/nFreqsPerIFS) * componentLossFactor;

				double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df, alpha);
				Complex[] fourierComp = getFourierComponentsForFrequencies(alpha, nuMin, nuMax, df); 
				// double[] nFrequencies = null;

				for (int i = 0; i < nFreqsPerIFS.length; i++) {
					fileLabel = "_si" + si + "_clf" + clf + "_nf" + nFreqsPerIFS[i];

					for (int j = 0; j < frequencies.length/((double) nFreqsPerIFS); j++) {						
						frequencies[j] = frequencies[nFreqsPerIFS[i]*j];
						fourierComp[j] = fourierComp[nFreqsPerIFS[i]*j];						
					}

					double[] rates = getRatesFromFourierComponents(fourierComp);
					double scalingFactor = mean/BasicStats.getMean(rates);

					for ( int i=0; i < rates.length; i++ ) {
					    rates[i] *= scalingFactor;
					}


					// Copied from RedNoiseGenerator
					double[] rateBins = new double[timmerRates.length];
				 	for (int i = 0; i < rateBins.length; i++) {
				 		rateBins[i] = i * duration/((double) timmerRates.length);
				 	}

				 // 	// Writes Time Domain Signal
				 // 	String fileName = "TDS-cr_" + meanRate + "d_" + duration + "si_" + alpha + "clf_" + componentLossFactor + "nf_" + nFreqsPerIFS + "b_" + nTimeBins;
				 // 	String[] header = AsciiDataFileWriter.makeHeader("Time-Domain Signal", ""); // labels axes
					// AsciiDataFileWriter writer = new AsciiDataFileWriter(fileName + ".qdp"); // Define filename
					// writer.writeData(header, rateBins, timmerRates);
					// System.out.println("\nWriting TimeDomainSignal as " + fileName + ".qdp\n");

					//  Compare values of minNBins and nTimeBins
					if ( minNBins != nTimeBins ) {
					    logger.warn("Closest power of 2 greater than specified number of time bins ("+nTimeBins+") not equal to number of bins required for effective Nyquist resolution ("+minNBins+")");
					}
				// 	if ( nTimeBins < minNBins ) {
				// 	    logger.warn(nTimeBins+" < "+minNBins+": This would cause loss of resolution: Resampling.");
				// 	    double[] oldBinEdges = BinningUtils.getBinEdges(0, duration, timmerRates.length);
				// 	    double[] newBinEdges = BinningUtils.getBinEdges(0, duration, nTimeBins);
				// 	    timmerRates = Resampler.resample(timmerRates, oldBinEdges, newBinEdges);	    
				// 	}

					//  Define the number of events
					Poisson poisson = new Poisson(0, engine); // defines a Poisson object
					int nevents = (new Double(meanRate*duration)).intValue();
					// TEST: remove randomness
					// nevents = poisson.nextInt(nevents); // takes a random value from the poisson distribution with mean nEvents
					// TEST

					//  Draw arrival times from the CDF of rates - see 'WeakPeriodicSignalsAndRedNoise' pg 42 for explanation
					double tzero = 0;
					double tkBinTime = duration/nTimeBins; // bin width
					Histogram1D lcHisto = Converter.array2histo("light curve", tzero, tkBinTime, timmerRates);
					Histogram1D cdfHisto = DistributionFunc.getCDFHisto(lcHisto);
					double[] times = DistributionFunc.getRandom(cdfHisto, nevents);
					Arrays.sort(times);


					// TEST: unnecessary statements
					// Adjust actual duration to specified duration
					// System.out.println("nEvents:\t" + nevents + "\ntimes.length\t" + times.length);
				 	// nevents=times.length; // ISSUE: this isn't necessary, they are already equal (see from uncommenting above line)
				 	// times[nevents-1] = times[0] + duration; // Prevents aioob on line 72 of PeriodogramGenerator (added catch block inside PeriodogramGenerator)
				 	// System.out.println("Final event time:\t" + times[nevents-1] );
				 	// TEST

					double actualMean = times.length/duration; // ISSUE: should be very close to actual mean given it was previously scaled to be the actual mean, so this is basically redundant.
					// TEST: how close is data mean to the input countRate? Answer: Very
					// System.out.println("dataMeanCountRate:\t" + actualMean);


					// Copied from Periodogram Generator

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
				}
				//

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						fileLabel = "_si" + si + "_clf" + clf + "_nf" + nFreqsPerIFS[i];
						
						// In Progress
						double nuMin = 1/(duration); // should be 128*duration
						double nuMax = mean;
						double df = (nuMin/nFreqsPerIFS) * componentLossFactor
		
						double[] frequencies = PeriodogramUtils.getFourierFrequencies(nuMin, nuMax, df, alpha);
				
						//

						periodogram = PeriodogramGenerator.periodogramGenerator(d, si, cr, b, nFreqsPerIFS[i], clf);
						runningAverageIndex += MinMaxFreq.maxFreq(periodogram, fileLabel);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;
					x[i] = (double) xInt[i]; // converts int to double to agree with other cases	
					System.out.println("nFreqsPerIFS: " + nFreqsPerIFS[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
				}
				y = maxFreqs;
				break;

			case "countRate-maxFreq": 
				x = countRate; 
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, countRate[i], b);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;
					System.out.println("countRate: " + countRate[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;				
				}
				y = maxFreqs;			
				break;

			case "duration-maxFreq": 
				x = duration; 
				maxFreqs = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(duration[i], si, cr, b);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;
					System.out.println("Duration: " + duration[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAverageIndex = 0;
				}
				y = maxFreqs;			
				break;

			case "duration-index": // plots estimatedIndex against duration
				x = duration;
				estimatedIndices = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(duration[i], si, cr, b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("Duration: " + duration[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
				}
				y = estimatedIndices;
				break;

			case "indices-index": // plots estimatedIndex against duration
				x = indices;
				estimatedIndices = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = ParameterEstimator.indexEstimator(d, indices[i], cr, b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("Indices: " + indices[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
				}
				y = estimatedIndices;
				break;

			case "countRate-index": // plots estimatedIndex against countRate
				x = countRate;
				estimatedIndices = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, si, countRate[i], b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("countRate: " + countRate[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
				}
				y = estimatedIndices;
				break;

			case "nBins-index": // plots estimatedIndex against nBins
				xInt = nBins; 
				x = new double[xInt.length];
				estimatedIndices = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, si, cr, nBins[i]);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;	
					x[i] = (double) xInt[i]; // converts int to double to agree with other cases
					System.out.println("nBins: " + nBins[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
				}	
				y = estimatedIndices;		
				break;

			// err cases plot on y |est-act|/act (est = estimtated index; act = actual index)
			case "duration-indexerr":
				x = duration;
				estimatedIndices = new double[x.length];
				estimatedIndexErrors = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(duration[i], si, cr, b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("Duration: " + duration[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
					estimatedIndexErrors[i] = Math.abs(estimatedIndices[i] - si)/si;
				}

				y = estimatedIndexErrors;
				break;

			case "indices-indexerr":
				x = indices;
				estimatedIndices = new double[x.length];
				estimatedIndexErrors = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, indices[i], cr, b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("Index " + indices[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
					estimatedIndexErrors[i] = Math.abs(estimatedIndices[i] - indices[i])/indices[i];
				}

				y = estimatedIndexErrors;
				break;

			case "countRate-indexerr":
				x = countRate;
				estimatedIndices = new double[x.length];
				estimatedIndexErrors = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, si, countRate[i], b);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("countRate: " + countRate[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
					estimatedIndexErrors[i] = Math.abs(estimatedIndices[i] - si)/si;
				}

				y = estimatedIndexErrors;
				break;

			case "nBins-indexerr":
				xInt = nBins; 
				x = new double[xInt.length];
				estimatedIndices = new double[x.length];
				estimatedIndexErrors = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, si, cr, nBins[i]);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					x[i] = (double) xInt[i]; // converts int to double to agree with other cases
					System.out.println("nBins: " + nBins[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
					estimatedIndexErrors[i] = Math.abs(estimatedIndices[i] - si)/si;
				}

				y = estimatedIndexErrors;
				break;

			case "nBins-maxPower": 
				xInt = nBins; 
				x = new double[xInt.length];
				maxPower = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						AveragePeriodogram avgPeriodogram = PeriodogramGenerator.periodogramGenerator(d, si, cr, nBins[i]);
						runningAveragePower += avgPeriodogram.getPowers()[0]; // returns the power value for the lowest frequency
					}
					maxPower[i] = runningAveragePower/totalCycles;	
					x[i] = (double) xInt[i]; // converts int to double to agree with other cases
					System.out.println("nBins: " + nBins[i] + "\t\tmaxPower: " + maxPower[i]);
					runningAveragePower = 0;
				}	
				y = maxPower;		
				break;

			case "lognBins-maxPower": 
				xInt = nBins; 
				x = new double[xInt.length];
				maxPower = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						AveragePeriodogram avgPeriodogram = PeriodogramGenerator.periodogramGenerator(d, si, cr, nBins[i]);
						runningAveragePower += avgPeriodogram.getPowers()[0]; // returns the power value for the lowest frequency
					}
					maxPower[i] = runningAveragePower/totalCycles;	
					x[i] = Math.log( (double) xInt[i]); // converts int to double to agree with other cases
					System.out.println("nBins: " + nBins[i] + "\t\tmaxPower: " + maxPower[i]);
					runningAveragePower = 0;
				}	
				y = maxPower;		
				break;

			case "logCountRate-logMaxFreq": 
				x = countRate;
				maxFreqs = new double[x.length];
				double[] logY = new double[x.length];
				double[] logX = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						runningAverageIndex += MinMaxFreq.maxFreq(d, si, countRate[i], b);
					}
					maxFreqs[i] = runningAverageIndex/totalCycles;
					logY[i] = Math.log( maxFreqs[i] );
					logX[i] = Math.log( x[i] ); 
					System.out.println("countRate: " + countRate[i] + "\t\tmaxFreq: " + maxFreqs[i]);
					runningAveragePower = 0;
				}
				x = logX;
				y = logY;

				System.out.println("\nLeast-squares fitting according to a power law...");
				fit = LeastSquaresFitter.leastSquaresFitLine(logX, logY);	
				System.out.println("\nGradient = " + fit[0] + "\nGradient error = " + fit[1] + "\nIntercept = " + fit[2] + "\nIntercept error = " + fit[3]);
	
				// System.out.println("\n\ndirect power law fitting...");
				// System.out.println("\nLeast-squares fitting according to a power law...");
				// fit = LeastSquaresFitter.leastSquaresFitPowerLaw(x, maxFreqs);	
				// System.out.println("\nIndex = " + fit[0] + "\nIndex error = " + fit[1] + "\nNormalisation = " + fit[2] + "\nNormalisation error = " + fit[3]);

				break;

			// edited changing type of correction factors from int to double
			case "correctionFactors-indexerr":
				x = correctionFactors; 
				estimatedIndices = new double[x.length];
				estimatedIndexErrors = new double[x.length];

				for (int i = 0; i < x.length; i++) {
					for (int cycles = 0; cycles < totalCycles; cycles++) {
						double estimatedIndex = IndexEstimator2.indexEstimator(d, si, cr, b, correctionFactors[i]);
						runningAverageIndex += estimatedIndex;
					}
					estimatedIndices[i] = runningAverageIndex/totalCycles;
					System.out.println("correctionFactors: " + correctionFactors[i] + "\t\tIndex estimate: " + estimatedIndices[i]);
					runningAverageIndex = 0;
					estimatedIndexErrors[i] = Math.abs(estimatedIndices[i] - si)/si;
				}

				y = estimatedIndexErrors;
				break;


			default: System.out.println("Not a valid plotting quantity"); System.exit(1);
		}

		// Writes plot as .qdp file
		String[] header = AsciiDataFileWriter.makeHeader(plot, ""); // labels axes

		// Decide on fileName
		// System.out.println("Enter QDP File Name: ");
		// String fileName = reader.next();

		// automatically generate fileName according to input
		String fileName = plot + "_d" + Double.toString(d) + "_si" + Double.toString(si) + "_cr" + Double.toString(cr) + "_b" + Integer.toString(b) + "_tc" + Double.toString(totalCycles) + "_cf" + Double.toString(cf) + "_nf" + Double.toString(nf);
		if (deterministic == 1) {
			fileName += "DET";
		}

		System.out.println(fileName);

		AsciiDataFileWriter writer = new AsciiDataFileWriter(fileName + ".qdp"); // Define filename

		//TEST
		//System.out.println(Arrays.toString(x) + "\n\n" + Arrays.toString(y));
		//TEST

		// if (x == null) {
		// 	writer.writeData(header, xInt, y); // Writes the data as a .qdp file. Note this method can take int[] or double[]	
		// }
		// else {
		// 	writer.writeData(header, x, y); 
		// }

		writer.writeData(header, x, y);

		System.out.println("To look at graph, use:\nqdp '" + fileName + "'.qdp\nthen type EXIT in command line to exit plot");


	}

}