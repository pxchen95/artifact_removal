# artifact_removal

## Explanation of common terms
* "segments" refer to segments of contiguous data on which stimulation is on and the stimulation parameters are **constant** (e.g., constant amplitude, fundamental frequency, shape, etc.).
* "gaps" refer to what is in between segments, whether that be missing data, data with stimulation off, data where the stimulation parameters do not match those of the stimulation in the segments, etc.

## Notation
*Included in comments in the code; typical notation with associated definitions.*
* fs  = sample rate in Hz
* n   = number of gaps in the data (or equivalently, total number of segments - 1)
* N_i = number of samples in i-th segment
