function [f_vec] = convert_cellarray_to_vector(f_cell, samp_shift, paddingtype)
% Convert signal stored as cell array to a vector with gaps padded using
% paddingtype
%
% INPUTS:
%      f_cell:      1 x numSegments cell array,
%                        f{i} = 1 x N_i, samples of segment i
%      samp_shift:  1 x numSegments-1 vector, 
%                        samp_shift(i) = length of gap in # of samples b/t segments i and i+1
%      paddingtype: 0 or NaN, fill gap using paddingtype
%
% OUTPUT:
%      f_vec: 1 x (sum(N_i) + sum(samp_shift)) vector, 
%                        samples of signal with gaps included

numSegments = length(f_cell);
f_vec = [];
for i = 1:numSegments
    f_vec = [f_vec f_cell{i}];
    if i < numSegments
        f_vec = [f_vec zeros(1,samp_shift(i))*paddingtype];
    end
end

end