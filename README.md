# artifact_removal

## Common Terms
* "segments" refer to segments of contiguous data on which stimulation is on and the stimulation parameters are **constant** (e.g., constant amplitude, fundamental frequency, shape, etc.).
* "gaps" refer to what is in between segments, whether that be missing data, data with stimulation off, data where the stimulation parameters do not match those of the stimulation in the segments, etc.

## Notation
*Included in comments in the code; typical notation with associated definitions.*
* fs  = sample rate in Hz
* n   = number of gaps in the data (or equivalently, total number of segments - 1)
* N_i = number of samples in i-th segment

## User Recommendations
* For better runtime, use at most ~1e4 total number of samples (across all segments). Preliminary tests showed that while using more than 1e4 samples usually resulted in more accurate frequency/phase shift estimates, the gain in accuracy did not significantly outweigh the loss in computational speed.
* If your data consists of only one segment, insert an artificial gap of 0 in the middle of the segment. Preliminary tests showed that estimating the frequency and a phase shift jointly is typically less sensitive to choice of initialization than searching for the frequency alone.
* Initialize newton_refinement_using_g.m with the frequency and phase shift estimates found using newton_rand_init.m. newton_refinement_using_g.m is very sensitive to initialization, but preliminary tests showed newton_rand_init.m usually provides a sufficiently good initialization.
* Note that segments are input into the functions as **cell arrays**.

## Brief Description of MATLAB Functions/Scripts
We group the functions/scripts as follows:
### Algorithm 1
*Estimates the frequency and phase shifts by maximizing the energy using Newton's ascent.*
* newton_rand_init.m: Runs Newton's ascent using uniform random initialization
* newton_ascent.m: Uses Newton's ascent to maximize the energy with respect to frequency and phase shifts
* remove_artifact.m: Uses simple harmonic regression to reconstruct and remove the artifact (optional)
* calc_E_multtshifts.m: Computes energy, its gradient, and its Hessian
* mod_cholesky.m: Implements modified Cholesky decomposition to ensure that a particular matrix in Newton's ascent is negative definite
* backtracking_linesearch: Computes the stepsize for Newton's ascent
 
### Algorithm 2
*Refines the frequency and phase shift estimates from Algorithm 1 by solving a least squares problem that minimizes over frequency, phase shift, and amplitudes through a combination of harmonic regression and Newton's descent.*
* newton_refinement_using_g.m: Runs Newton's descent to solve the least squares problem to refine the estimates and reconstruct/remove the artifact
* remove_artifact_ver_g.m: Computes the objective function in the least squares problem, its gradient, and its Hessian 
* backtracking_linesearch_for_g.m: Computes the stepsize for Newton's descent (currently unused)

### Other
* demo/demo.m: Demonstration of how to use the included functions with examples of different ways to visualize the results
* demo/1segment_example.mat: 1 long segment of data, where artifact is sum of sinusoids with 5 harmonics, true underlying signal is Gaussian noise
* demo/multsegment_example.mat: 10 short segments of data, where artifact is sum of sinusoids with 5 harmonics, true underlying signal is Gaussian noise
* convert_cell_array_to_vector.m: Converts cell array to vector with gaps specified by input *samp_shift*
* convert_vector_to_cellarray.m: Converts vector to cell array with lengths of segments specified by input *N* and lengths of gaps specified by input *samp_shfit*
