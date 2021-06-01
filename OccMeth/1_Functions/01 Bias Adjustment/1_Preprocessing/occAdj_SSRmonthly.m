function [outho, ouths, outfs, th] = occAdj_SSRmonthly(xho, xhs, xfs)
%   OCCADJ_SRR This function implements a precipitation occurrence bias adjustment 
%   method using the stochastic singularity removal technique
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   This function implements the stochastic singularity removal as
%   proposed by Vrac, Noël and Vautard (2016) in "Bias correction of
%   precipitation through Singularity Stochastic Removal: Because
%   occurences matter"
%
%   Inputs:
%       xfs: the future simulations, to be corrected, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xho: the historical observations, idem
%       xhs: the historical simulations, idem
%   Outputs:
%       outho: matrix of the changed historical observations
%       ouths: matrix of the changed historical simulations
%       outfs: matrix of the changed future simulations
%       th: threshold calculated in this function
%
%   Last update by J. Van de Velde on 08/12/'20

%% Implementation

% Setup

outfs = xfs(:, 4);
ouths = xhs(:, 4);
outho = xho(:, 4);

th = zeros(1,12);

% Monthly loop

for m = 1:12 %For each month
    idm = find(xfs(:,2) == m); %Gets all the indices of the month m
    
    %Selection of the rows belonging to the month m
    xfsm = xfs(xfs(:,2) == m ,:);
    xhsm = xhs(xhs(:,2) == m,:);
    xhom = xho(xho(:,2) == m,:);
    
    ndays = length(xhom);
    
    % Determine th
    % th is determined by the smallest wet day precipitation depth
    
    thm = min([min(xhom(xhom(:,4)>0, 4)), min(xhsm(xhsm(:, 4)>0, 4)), min(xfsm(xfsm(:, 4)>0, 4))]);
    
    % Loop over the days
    
    for i= 1:ndays
        if xhom(i, 4) < thm
            randomnmbr = rand(1);
            xhom(i,4) = randomnmbr*thm; %Random selection of a day out of the interval between 0 and the threshold
        end
        if xhsm(i, 4) < thm
            randomnmbr = rand(1);
            xhsm(i, 4) = randomnmbr*thm;
        end
        if xfsm(i,4) < thm
            randomnmbr = rand(1);
            xfsm(i,4) = randomnmbr*thm;
        end
    end
    
    % Arrange in original order in total
    for i = 1:length(idm)
        outho(idm(i)) = xhom(i,end);
        ouths(idm(i)) = xhsm(i,end);
        outfs(idm(i)) = xfsm(i,end);
    end
    
    th(m) = thm;
    
end
end

