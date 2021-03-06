function [out_SSR] = postprocessingSSR(out, pth, p)
%   postprocessingSSR This function post-processes the quantile mapping output
%   in case the SSR wet day frequency BA method was used.
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This file postprocesses the bias adjusted data using the threshold
%   implemented in SSR. All values of the bias adjusted data that are
%   below this threshold, are changed to zero and the statistics are
%   calculated.
%
%   Inputs:
%       out: intensity-adjusted time series
%       pth: threshold calculated in the first correction step by SSR
%       p: column number of precipitation in 'out'
%   Output:
%       out_SSR: fully adjusted time series
%
%   Last update by J. Van de Velde on 10/12/'20

%% Set-up

n = size(out,2);
ntime = length(out{1});
out_SSR = out;

%% Loop

for i = 1:n
    for j=1:ntime
        m = out{i}(j,2);
        if out{i}(j,p) < pth(m)
            out_SSR{i}(j,p) = 0;
        end
    end

end

