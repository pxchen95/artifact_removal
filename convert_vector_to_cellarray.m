function [f_cell] = convert_vector_to_cellarray(f_vec, N, samp_shift)
% Convert signal stored as vector with gaps padded with NaNs to cell array
%
% INPUTS:
%      f_vec:       1 x (sum(N_i) + sum(samp_shift)) vector, 
%                        samples of signal with gaps included
%      N:           1 x numSegments vector, 
%                        N(i) = N_i, sample length of segment i
%      samp_shift:  1 x numSegments-1 vector, 
%                        samp_shift(i) = length of gap in # of samples b/t segments i and i+1
%
% OUTPUT:
%      f_cell: 1 x numSegments cell array, f{i} = 1 x N_i, samples of segment i

f_cell = {};
numSegments = length(N);
ind1 = 1;
for i = 1:numSegments
    ind2 = ind1 + N(i) - 1;
    f_cell{i} = f_vec(ind1:ind2);
    if i < numSegments
        ind1 = ind2 + 1 + samp_shift(i);
    end
end

end