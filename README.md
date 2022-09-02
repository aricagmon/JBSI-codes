# JBSI-codes
These are MathCad and MatLab codes to compute the JBSI (Jitter-Based Synchrony Index) of a simultaneously recorded pair of spike trains.
See the original publication (Agmon A (2012) A novel, jitter-based method for detecting and measuring spike synchrony and quantifying temporal firing precision. Neural systems & circuits 2:5. doi:10.1186/2042-1001-2-5) for details on the algorithm.

There are two MathCad codes, one to calculate a single JBSI value based on given synchrony window and shift values, the other to calculate a matrix of JBSI  for a range of synchrony windows and shift values. These routines can be tried with any of the 3 Excel data files, although the matrix calculation may take a long time when used with the longer data files. 

The MatLab function was written by Craig Atencio, and is essentially a translation of the single-JBSI MathCad code. To try out this code, load the matlab data file. There will be two spike trains labeled a and b (ignore the variable labeled 'data'). Then issue the command:
>>syncData = jbsi_synchrony(a, b);
This will produce a list of synchrony and jitter windows that are being used, and a figure will appear containing the cross-correlogram and the jbsi for different windows.
syncData is a struct that holds the results. It holds all the calculations.
