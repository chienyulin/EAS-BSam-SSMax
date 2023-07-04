function [data, M, m] = scale_func(data,M,m)
%
% This function rescale the input data between -1 and 1
%
% INPUTS:
%   - data: the data to rescale
%   - M:    the maximum value of the output data
%   - m:    the minimum value of the output data
% OUTPUTS:
%   - data: the rescaled data
%
[Nb_s, Nb_b]=size(data);
if nargin==1
    M=max(data,[],1);
    m=min(data,[],1);
end

data = 2*(data-repmat(m,Nb_s,1))./(repmat(M-m,Nb_s,1))-1;
data(isnan(data))=0;
data(isinf(data))=0;
