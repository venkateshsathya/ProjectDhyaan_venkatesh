
function saveConnDataWithGrangerGoodChannelsDFTParametric(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange, connDataFolder, add_params)
% connMethod: 'granger_parametric' (pairwise AR), 'granger_bivariate',
% 'granger_multivariate', or an existing
% FieldTrip method such as 'ppc' or 'coh'. 'granger' uses multivariate GC.
% Select good electrodes before AR fitting, or after Fourier analysis for nonparametric GC.
% Saved matrices keep the full channel order, with bad rows/columns set to NaN.
% Subtract fitted 50/100 Hz sinusoids per trial/window before Fourier analysis.
% The connMethod suffix distinguishes parametric results; save add_params with results.

%Additional input parameters: add_params
% Parametric settings: mvar_order, mvar_method (bsmart/biosig),
% mvar_algorithm (BioSig only), mvar_demean (yes/no).
% Fits a separate bivariate model for every pair, pooled across trials.

if ~exist('stRange','var');         stRange = [0.25 1.25];              end
% % foilim_max = data.fsample/2;
% add_params.remove_bad_eletrod_flag = true;

numProtocols = length(protocolNameList);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder,add_params);
end
end

function saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder,add_params)

%connDataFolder = 'savedData'; % local project directory
%makeDirectory(connDataFolder);


[success, message] = mkdir(connDataFolder);
if ~success
    error('Could not create folder: %s', message);
end

disp(connDataFolder);
connDataFolder = fullfile(connDataFolder,subjectName);
makeDirectory(connDataFolder);

% Analysis file - this is where the final analyzed data is saved
analysisDetailsFileConn = fullfile(connDataFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']);

% ftData - data in fieldtrip format is saved in ftData
ftDataFileName = fullfile(ftDataFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
tmpData = load(ftDataFileName);
numGoodTrials = tmpData.numGoodTrials;

if numGoodTrials>0
    data = tmpData.data;

    cfg        = [];
    cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0]; % just to have data of same length in Pre, Post


    dataForFT_NoTimevals = data;
    % if isfield(dataForFT_NoTimevals, 'timeVals')
    %     dataForFT_NoTimevals = rmfield(dataForFT_NoTimevals, 'timeVals');
    % end

    dataPre = ft_redefinetrial(cfg, dataForFT_NoTimevals);

    %dataPre   = ft_redefinetrial(cfg,data);

    cfg        = [];
    cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
    dataPost   = ft_redefinetrial(cfg,dataForFT_NoTimevals);

    [connPre,freqPre] = getConnIndividualSubject(dataPre,connMethod,add_params);
    [connPost,freqPost] = getConnIndividualSubject(dataPost,connMethod,add_params);
else
    connPre=[]; connPost=[]; freqPre=[]; freqPost=[];
end

save(analysisDetailsFileConn,'connPre','connPost','freqPre','freqPost', ...
    'numGoodTrials','ftDataFileName','add_params');

end

function [conn,freqVals] = getConnIndividualSubject(data,connMethod,add_params)

[temp,fieldtripMethod]=getConnThisCondition(data,connMethod,add_params);
if any(strcmp(connMethod,{'granger_bivariate','granger_parametric'}))
    % Granger bivariate provides the output in a fashion of NC2*2 -
    % pairwise selected Granger FC and into 2 for two directions.
    % We need to move data to a format of 64 by 64 by 2. multivariate/PPC
    % etc by default provides FC ouput in this dimension. Downstream data
    % needs data in this dimension.
    % Extracts the pair of labels from temp.labelcmb, first column of pairLabels is
    % source, second is target.
    pairLabels = regexprep(temp.labelcmb,'\[.*$','');
    %Converts the labels in text such as Oz etc to numbers.
    [~,sourceIdx] = ismember(pairLabels(:,1),data.label);
    [~,targetIdx] = ismember(pairLabels(:,2),data.label);
    conn = zeros(numel(data.label),numel(data.label),numel(temp.freq));
    for k=1:size(pairLabels,1)
        conn(sourceIdx(k),targetIdx(k),:) = reshape(abs(temp.grangerspctrm(k,:)),1,1,[]);
    end
elseif strcmp(fieldtripMethod,'granger')
    % Restore the good-channel matrix to the original full channel order.
    [~,channelIdx] = ismember(temp.label,data.label);
    conn = nan(numel(data.label),numel(data.label),numel(temp.freq));
    conn(channelIdx,channelIdx,:) = abs(temp.grangerspctrm);
else
    conn = abs(temp.([fieldtripMethod 'spctrm']));
end
badElecs = data.badElecs;

conn(badElecs,:,:) = NaN;
conn(:,badElecs,:) = NaN;

freqVals = temp.freq;
end
function [conn_stat,fieldtripMethod] = getConnThisCondition(data,method,add_params)
% Clean each channel independently; keep the original trial/time metadata.
if add_params.linenoise_flag
    for iTrial=1:numel(data.trial)
        data.trial{iTrial} = ft_preproc_dftfilter(data.trial{iTrial}, ...
            data.fsample,add_params.line_noise_freq,'dftreplace',add_params.ln_remove_method);
    end
end
switch method
    case 'granger_parametric'
        % Fit each electrode pair directly in time, before evaluating spectra.
        if add_params.remove_bad_eletrod_flag
            goodChannels = true(numel(data.label),1);
            goodChannels(data.badElecs) = false;
            cfg = [];
            cfg.channel = data.label(goodChannels);
            data = ft_selectdata(cfg, data);
        end
        if numel(data.label)<2
            error('Parametric Granger requires at least two selected electrodes.');
        end
        validateattributes(add_params.mvar_order, {'numeric'}, ...
            {'scalar','integer','positive','finite'}, mfilename, 'add_params.mvar_order');
        if any(cellfun(@(x) size(x,2), data.trial)<=add_params.mvar_order)
            error('AR order must be smaller than the number of samples in every trial.');
        end
        validateattributes(add_params.foilim_max, {'numeric'}, ...
            {'scalar','positive','finite','<=',data.fsample/2}, mfilename, 'add_params.foilim_max');
        cfg = [];
        cfg.method = add_params.mvar_method;
        cfg.order = add_params.mvar_order;
        if strcmp(cfg.method,'biosig')
            cfg.mvarmethod = add_params.mvar_algorithm;
        end
        cfg.channelcmb = {'all','all'}; % Separate bivariate fits, not one joint model.
        cfg.keeptrials = 'no';          % Pool trials while respecting trial boundaries.
        cfg.demean = add_params.mvar_demean;
        cfg.zscore = 'no';
        cfg.taper = 'rectwin';
        mdata = ft_mvaranalysis(cfg, data);
        cfg = [];
        cfg.foi = 0:1:add_params.foilim_max;
        % Call the MVAR evaluator directly: the generic wrapper in this
        % FieldTrip version expects label, but pairwise models use labelcmb.
        freqPost = ft_freqanalysis_mvar(cfg, mdata);
    otherwise
        cfg              = [];
        cfg.method       = 'mtmfft';
        cfg.taper        = 'dpss';
        cfg.tapsmofrq    = 1;
        if(strcmp(method,'coh'))
            cfg.output = 'powandcsd';
        else
            cfg.output = 'fourier';
            cfg.keeptrials = 'yes';
        end

        cfg.foilim       = [0 add_params.foilim_max]; % Full spectrum for Granger factorization.
        freqPost = ft_freqanalysis(cfg, data);
        % Remove the data corresponding to bad electrodes
        if any(strcmp(method,{'granger_bivariate','granger_multivariate','granger'})) && add_params.remove_bad_eletrod_flag
            % Find the bad electrodes
            goodChannels = true(numel(data.label),1);
            goodChannels(data.badElecs) = false;
            cfg = [];
            %Assign the config param channel to those with letter names
            cfg.channel = data.label(goodChannels);
            % FreqPost has now data and labels trimmed to those of only good
            % electrodes
            freqPost = ft_selectdata(cfg, freqPost);
        end
end

if(strcmp(method,'coh'))
    TfreqPost = ft_checkdata(freqPost,'cmbrepresentation', 'full');
else
    TfreqPost = freqPost;
end

cfg=[];
switch method
    case {'granger_bivariate','granger_parametric'}
        cfg.method = 'granger';
        cfg.granger.sfmethod = 'bivariate';
        cfg.granger.conditional = 'no';
        cfg.channelcmb = {'all','all'};
    case {'granger_multivariate','granger'}
        cfg.method = 'granger';
        cfg.granger.sfmethod = 'multivariate';
        cfg.granger.conditional = 'no';
    otherwise
        cfg.method = method;
end
fieldtripMethod = cfg.method;
conn_stat = ft_connectivityanalysis(cfg, TfreqPost);
end
