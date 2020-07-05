# SimpleSVM
 SVMs are used here to determine whether the momentary actvity pattern of an ensemble of neurons can discriminate the passage of time within a 300 msec interval defined by 2 stimuli. SVM categorization was used specifically to test whether the ensemble of neurons showed activity patterns that could uniquely and reliably predict the correct time bin during the stimulus-free interval. 
 
 To Run:
 Load one dataset and then type 'SVMsimple_v8' (no quotes) in command window to run.
 Most interesting to run '3bin' data first, and then note the different result for the '6bin' dataset.
 
 Skills: 
 Demonstrates implentation of Multi-Class SVM, and more importantly use of single SVMs to provide deep insight regarding Multi-Class SVM result. Innovative use of bootstrapping approach is implented to control for data set partition issues depending on random assignment of samples to training or testing datasets, providing confidence in the results and interpretation.
 
 NOTE: Code to pull the data from the original files is included, but will not run in the absence of the original data files. Pre-processed structured sample data is included with this package. If loaded first, the code will skip the data pull and structuring section. The code is included here as an example, and was used to create the sample data sets.
 
 General Information:
 The passage of interval time is defined as discrete time bins, the size of which (50 or 100 ms) yields information regarding the resolution of the "neural clock". These interval time bins are the labels used for the SVM
 
 400 interval samples of spike data are used for each neuron. The sample intervals are split into training and test sets for the SVM, and are inherently fully balanced for each interval time bin. Data for each interval for each neuron was rescaled using min-max (0-1) transformation.
 
 In addition to a standard Multi-Class SVM built-in Matlab function as a global assessment of performance, I also used single SVMs for cross-bin comparisons to more closely examine which time bins were robustly encoded by neurons, and which time bins showed a decreased ability to discriminate interval time. Because the partiular random split of the data into training and test sets can probabistically lead to an unrepresentative result, the bin-x-bin comparisons were each bootstrapped to randomly resample the data into training and test sets and the SVM re-trained and re-tested iteratively to create a distribution reflecting the most accurate assessment of SVM performance independent of data assignment.
 
 Results and Interpretation:
 For 3 interval bins (100 msec each), the multi-class SVM suggests good performance with a low error rate. The single SVM analysis shows that the neural data can be used to decode interval time for all interval time bins with a high level of accuracy, well above chance levels. However, when the temporal resolution is increased to 50 msec (6 time bins for the 300 msec interval), the multi-class SVM indicates decreased performance and a higher level of information loss. The single SVM analysis reveals that as interval time proceeds, later time bins are not discriminated locally above chance levels, while early interval time bins are still discriminated well. The decrease in SVM performance with greater temporal resolution is specifically due to a decreased ability to discrimate momentary time late in the interval, and not poorer performance over all time bins.
 
 
 
 
