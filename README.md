# artifact_removal

Patents related to these algorithms have been provisionally filed.

Link to Github code: https://github.com/pxchen95/artifact_removal

## Common Terms
* "segments" refer to segments of contiguous data on which stimulation is on and the stimulation parameters are **constant** (e.g., constant amplitude, fundamental frequency, shape, etc.).
* "gaps" refer to what is in between segments, whether that be missing data, data with stimulation off, data where the stimulation parameters do not match those of the stimulation in the segments, etc. 
* "stimulation" and "artifact" are used interchangeably.
* "phase shift i" refers to the phase shift between the periodic artifact in the 0-th and i-th segment of data and corresponds to the i-th gap in the data.

## Notation
*Included in comments in the code; typical notation with associated definitions.*
* fs  = sample rate in Hz
* n   = number of gaps in the data (or equivalently, total number of segments - 1)
* N_i = number of samples in i-th segment

## User Recommendations
* For better runtime, use at most ~1e4 total number of samples (across all segments). Preliminary tests showed that while using more than 1e4 samples usually resulted in more accurate frequency/phase shift estimates, the gain in accuracy did not significantly outweigh the loss in computational speed. 
* If your data consists of only one segment, insert an artificial gap of 0 in the middle of the segment. Preliminary tests showed that estimating the frequency and a phase shift jointly is typically less sensitive to the choice of initialization than estimating the frequency alone.
* Initialize newton_refinement_using_g.m with the frequency and phase shift estimates found using newton_rand_init.m. newton_refinement_using_g.m is sensitive to initialization, but preliminary tests showed that newton_rand_init.m usually provides a sufficiently good initialization.
* Initialize newton_rand_init.m with the device setting of the frequency. Preliminary tests showed that the device setting of the frequency usually provides a sufficiently good initialization.
* Preliminary results show that Algorithm 1 (newton_refinement_using_g.m) provides more accurate frequency and phase shift estimates than Algorithm 2 (newton_rand_init.m), but we also include with Algorithm 2 some code implementing simple harmonic regression (remove_artifact.m), which can reconstruct/remove the artifact using the frequency and phase shift estimates from Algorithm 2 (as opposed to those from Algorithm 1). It was observed that in some cases (e.g., when the amplitude of the artifact is much higher than that of the underlying neural signal), the inaccuracies in the Algorithm 2 estimates led to "edge effects" in the reconstructed signal when using simple harmonic regression (i.e., larger errors near the ends of the segments of the reconstructed signal). However, in some cases, using Algorithm 2 + simple harmonic regression can have faster runtimes than (Algorithm 2 +) Algorithm 1.
* Note that segments are input into the functions as **cell arrays**.

## Brief Description of MATLAB Functions/Scripts
We group the functions/scripts as follows: 
### Algorithm 1
*Estimates the frequency and phase shifts by solving a least squares problem that minimizes over frequency, phase shifts, and amplitudes by jointly applying harmonic regression and Newton's descent. Corresponds to Algorithm 1 (artifact removal algorithm) in [1].*
* newton_refinement_using_g.m: Runs Newton's descent to solve the least squares problem and to refine the estimates and reconstruct/remove the artifact
* remove_artifact_ver_g.m: Computes the objective function in the least squares problem, its gradient, and its Hessian 
* backtracking_linesearch_for_g.m: Computes the stepsize for Newton's descent (currently unused)

### Algorithm 2 (Initialization Algorithm for Algorithm 1)
*Estimates the frequency and phase shifts by maximizing the energy using Newton's ascent. Corresponds to Algorithm 2 (initialization algorithm) in [1].* Also includes an option to use simple harmonic regression with the frequency and phase shift estimates from Algorithm 2 to reconstruct/remove the artifact.
* newton_rand_init.m: Runs Newton's ascent using uniform random initialization
* newton_ascent.m: Uses Newton's ascent to maximize the energy with respect to frequency and phase shifts
* remove_artifact.m: Uses simple harmonic regression to reconstruct and remove the artifact (optional)
* calc_E_multtshifts.m: Computes energy, its gradient, and its Hessian
* mod_cholesky.m: Implements modified Cholesky decomposition to ensure that a particular matrix in Newton's ascent is negative definite
* backtracking_linesearch: Computes the stepsize for Newton's ascent

### Other
* demo/demo.m: Demonstration of how to use the included functions with examples of different ways to visualize the results
* demo/1segment_example.mat: 1 long segment of data, where artifact is sum of sinusoids with 5 harmonics, true underlying signal is independent and identically distributed (iid) Gaussian noise
* demo/multsegment_example.mat: 10 short segments of data, where artifact is sum of sinusoids with 5 harmonics, true underlying signal is iid Gaussian noise
* paper_examples/run_paper_examples.m: Runs numerical examples 1-3 in [1]
* paper_examples/example1_SingleSegmentArtifactOnly.mat: 1 segment of data, where artifact is sum of sinusoids with 5 harmonics, no underlying signal
* paper_examples/example2_SingleSegmentChirp.mat: 1 segment of data, where artifact is sum of sinusoids with 5 harmonics, true underlying signal is a chirp
* paper_examples/example3_ManySegmentsAliased.mat: 10 segments of data, where artifact is sum of sinusoids with 5 harmonics and the fundamental frequency is aliased due to low sampling rate, true underlying signal is a simulated neural signal (computed as the sum of short snippets of sinusoidal waves with random frequencies/random lengths and demeaned) plus iid Gaussian noise
* convert_cell_array_to_vector.m: Converts cell array to vector with gaps specified by input *samp_shift*
* convert_vector_to_cellarray.m: Converts vector to cell array with lengths of segments specified by input *N* and lengths of gaps specified by input *samp_shift*

## References
[1] P. Chen, et al, "Estimation of Periodic Signals with Applications to Deep Brain Stimulation," *preprint*, 2022.

Note: This paper is available on this Github as EstimationPeriodicSignalDBS.pdf.
