function BiasAdjustment(xho, xhs, xfs, var, methodocc, methodint, n, modeltype)
%   BiasAdjustment This function implements the bias adjustment according
%   to the input given.
% 
%   This function is launched in the b_configurationBiasAdjustment config
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   Inputs: 
%       xho: vector that contains the original observations for the current/control period
%       xhs: vector that contains the original GCM output for the current/control period 
%       xfs: vector that contains the original GCM output for the future period
%       var: specifies the input variable: "P" for Precipitation, "E" for Evapoporation or "T" for Temperature
%       methodocc: 'none', 'tda', 'ssr'
%       methodint: 'qdm' 
%       n: ensemble members or repetitions base number
%       modeltype: string with the name of the model run used
%
%   Last update by J. Van de Velde on 10/12/'20



%% Set-up

% Adding save location
save_bias = strcat('E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\',modeltype,'_', methodocc,'_',methodint); %Bias correction path, file named is created by combining the specifications

% Data preparation
ndays = length(xho);
timef = xfs(:,1:3);
timeh = xho(:,1:3);

% Making selections of the variables
% The datasets are structured as follows: Y M D VAR VAR VAR, so to select
% the variable, its location in 'var' + 3 is needed.
% This is based on the structure of the Uccle dataset, and might need
% changes if another dataset is used.

for i = 1:length(var)
    if strcmp(var{i},'T') == 1
        t = i+3;
        tho = xho(:,[1:3,t]); %All rows, columns 1:3 (Y M D) + the variable 
        ths = xhs(:,[1:3,t]);
        tfs = xfs(:,[1:3,t]);
    elseif strcmp(var{i},'E') == 1
        e = i+3;
        eho = xho(:,[1:3,e]);
        ehs = xhs(:,[1:3,e]);
        efs = xfs(:,[1:3,e]);
    else
        p = i+3;
        pho = xho(:,[1:3,p]);
        phs = xhs(:,[1:3,p]);
        pfs = xfs(:,[1:3,p]);
    end
end

%% Preprocessing: Occurrence bias adjustment

% In this part, occurrence bias adjustment takes place

if strcmp(methodocc, 'none') == 0
    switch methodocc
        case 'tda'
            %Initialize dry-wet corrected vectors
            pfsdw = nan(ndays,n+3); %20 repetitions, extra column for each member
            pfsdw(:,1:3) = timef;
            phsdw = nan(ndays,n+3);
            phsdw(:,1:3) = timeh;
            phodw = nan(ndays, n+3);
            phodw(:, 1:3) = timeh;
            
            for i= 1:n
                %Adjustment
                [pfsdw(:,i+3), phsdw(:,i+3)] = occAdj_TDA(pfs,pho,phs);
                %Repetition of unchanged variables
                phodw(:,i+3) = pho(:,4); % Unchanged, so repetition of the same column
            end
            
        case 'ssr'
            %Initialize dry-wet corrected vectors
            pfsdw = nan(ndays,n+3); %20 repetitions, extra column for each member
            pfsdw(:,1:3) = timef;
            phsdw = nan(ndays,n+3); %20 repetitions, extra column for each member
            phsdw(:,1:3) = timeh;
            phodw = nan(ndays, n+3); %20 repetitions, extra column for each member
            phodw(:, 1:3) = timeh;
            pth = nan(12, n);
            
            for i= 1:n
                [phodw(:,i+3), phsdw(:,i+3), pfsdw(:,i+3), pth(:,i)] = occAdj_SSRmonthly(pho,phs,pfs);
                %[phodw2(:,i+3), phsdw2(:,i+3), pfsdw2(:,i+3), pth2] =
                %occAdj_SSRor(pho,phs,pfs); Changed SSR to work monthly, tested
                %whether this had an influence. There is an influence, but it
                %is very small.
            end
    end
else
    %'No explicit adjustment' of precipitation
    pfsdw = [timef pfs(:, 4)];
    phsdw = [timeh phs(:, 4)];
    phodw = [timeh pho(:, 4)];
end

%% Intensity Bias adjustment

% In this part, the intensity bias adjustment takes place

%Initialize rescaled vectors
pfc = nan(ndays,n+3);
pfc(:,1:3) = timef;

switch methodint
    case 'qdm'
        n2 = (size(phodw,2)-3); %Repetition number based on occurrence bias adjustment
        for i = 1:n2 % For each of the repetitions made
        [~, pfc(:,i+3)] = QDM(phodw(:,[1:3,i+3]),phsdw(:,[1:3,i+3]),pfsdw(:,[1:3,i+3]), 2);
        end
        [~, efc] = QDM(eho,ehs,efs, 2);
        [~, tfc] = QDM(tho,ths,tfs, 1); %QDM and occurrence bias adjustment aren't needed for temperature data
        
        %Saving
        out = cell(1,n2); %Cell size depending on repetitions needed
        for i = 1:n2 %Save everything, because some variables will be overwritten
            out{i} = nan(size(xho));
            out{i}(:,1:3) = timef;
            out{i}(:,t) = tfc(:,end);
            out{i}(:,p) = pfc(:,i+3);
            out{i}(:,e) = efc(:,end);
        end
        save(strcat(save_bias, '_results'),'out')
    
end

%% Post-processing

% Post-processing steps, like in SSR, end up here
switch methodocc
    case 'ssr'
         [out] = postprocessingSSR(out, pth, p);
        save(strcat(save_bias, '_results'),'out') %Overwrites previous version, though it might be interesting to compare before and after postprocessing
end

end
