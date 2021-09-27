function [M, tau] = mod_cholesky(A)
% Returns matrix M = A - tau*I, where I is the identity matrix, such that M
% is a negative definite matrix
%
% INPUTS:
%    A: n x n matrix
%
% OUTPUTS: 
%    M:   n x n matrix
%    tau: scalar

normA = norm(A, 'fro'); % froebenius norm
n = size(A,1);

if normA <= 1e-6
    normA = 1e-6;
end

if max(diag(A)) < 0 % A is already negative definite
    tau = 0;
else
    tau = normA;
end

maxEigValue = max(eig(A));
maxEigValuePlusTau = maxEigValue - tau;

% find a tau that makes M = A - tau*I negative definite
while maxEigValuePlusTau >= 0 
    tau = max(2*tau, 0.5*normA);
    maxEigValuePlusTau = maxEigValue - tau;
end

M = A - tau*eye(n);

end