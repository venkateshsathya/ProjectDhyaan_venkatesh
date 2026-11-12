function saveConnDataWithGrangerGoodChannelsDFTParametricFast(subjectName, protocolNameList, badEyeCondition, badTrialVersion, ftDataFolder, connMethod, stRange, connDataFolder, add_params)
% Input options are set in runSaveFTAndConnDataParametricFast.m.
% Granger: parametric/nonparametric, with bivariate/multivariate models.
% Other methods: coherence, PPC and other Fourier-based FieldTrip methods.
% keeptrials = 'no' is implemented; 'yes' is a placeholder for future work.
% Saved matrices retain the original channel order and bad rows/columns are NaN.

if ~exist('stRange', 'var')
    stRange = [0.25 1.25];               % Use the existing default Post window.
end

% Defaults for the new options allow existing scripts to call this function.
defaultParams = struct('granger_method', 'nonparametric', ...
    'granger_model', 'multivariate', 'keeptrials', 'no', ...
    'taper', 'dpss', 'tapsmofrq', 1, 'mvar_freqstep', 1, ...
    'mvar_zscore', 'no', 'mvar_taper', 'rectwin');
optionNames = fieldnames(defaultParams);
for iOption = 1:numel(optionNames)
    optionName = optionNames{iOption};
    if ~isfield(add_params, optionName)
        add_params.(optionName) = defaultParams.(optionName);
    end
end

% Map the explicit Granger names to their calculation options.
analysisMethod = connMethod;
switch connMethod
    case 'granger_parametric'
        analysisMethod = 'granger';
        add_params.granger_method = 'parametric';
        add_params.granger_model = 'bivariate';
    case 'granger_parametric_multivariate'
        analysisMethod = 'granger';
        add_params.granger_method = 'parametric';
        add_params.granger_model = 'multivariate';
    case 'granger_nonparametric_bivariate'
        analysisMethod = 'granger';
        add_params.granger_method = 'nonparametric';
        add_params.granger_model = 'bivariate';
    case 'granger_nonparametric_multivariate'
        analysisMethod = 'granger';
        add_params.granger_method = 'nonparametric';
        add_params.granger_model = 'multivariate';
end

% This switch applies to every method and every Granger model below.
switch add_params.keeptrials
    case 'no'
        % Pool the retained trials, as in the existing pooled calculation.
    case 'yes'
        % Placeholder: add a separate FC calculation for each trial here later.
        error('keeptrials = yes is a placeholder. Use add_params.keeptrials = ''no'' for now.');
end

% Compute and save each protocol separately for this subject.
for iProtocol = 1:numel(protocolNameList)
    protocolName = protocolNameList{iProtocol};
    saveConnDataSingleProtocol(subjectName, protocolName, badEyeCondition, ...
        badTrialVersion, ftDataFolder, connMethod, analysisMethod, stRange, ...
        connDataFolder, add_params);
end
end

function saveConnDataSingleProtocol(subjectName, protocolName, badEyeCondition, badTrialVersion, ftDataFolder, connMethod, analysisMethod, stRange, connDataFolder, add_params)

% Keep the subject folders and include the selected method name in each filename.
mkdir(connDataFolder);
disp(connDataFolder);
connDataFolder = fullfile(connDataFolder, subjectName);
makeDirectory(connDataFolder);
analysisDetailsFileConn = fullfile(connDataFolder, ...
    [protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']);

% Load the FieldTrip data and the retained trial count.
ftDataFileName = fullfile(ftDataFolder, subjectName, ...
    [protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
tmpData = load(ftDataFileName);
numGoodTrials = tmpData.numGoodTrials;

if numGoodTrials > 0
    data = tmpData.data;

    % Optional resampling before selecting the Pre and Post windows.
    if isfield(add_params, 'resamplefs') && ~isempty(add_params.resamplefs) ...
            && add_params.resamplefs < data.fsample
        badElecs = data.badElecs;
        cfg = [];
        cfg.resamplefs = add_params.resamplefs;
        cfg.method = 'resample';         % Includes anti-alias filtering.
        cfg.detrend = 'no';
        data = ft_resampledata(cfg, data);
        data.badElecs = badElecs;        % Retain the bad-electrode indices.
    end
    
    % Pre and Post contain the same number of samples, as before.
    cfg = [];
    cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0];
    dataPre = ft_redefinetrial(cfg, data);
    cfg = [];
    cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
    dataPost = ft_redefinetrial(cfg, data);

    % Compute FC separately for the Pre and Post windows.
    [connPre, freqPre] = getConnIndividualSubject(dataPre, analysisMethod, add_params);
    [connPost, freqPost] = getConnIndividualSubject(dataPost, analysisMethod, add_params);
else
    % Retain the existing empty outputs when there are no good trials.
    connPre = []; connPost = []; freqPre = []; freqPost = [];
end

% For Granger, rows are sources and columns are targets; frequency is dimension 3.
% Coherence and PPC are symmetric, so the two channel directions are equivalent.
connDimord = 'source_target_freq';
save(analysisDetailsFileConn, 'connPre', 'connPost', 'freqPre', 'freqPost', ...
    'numGoodTrials', 'ftDataFileName', 'add_params', 'connDimord', '-v7.3');
end

function [conn, freqVals] = getConnIndividualSubject(data, method, add_params)

% Retain the original labels before selecting good channels.
allLabels = data.label;
badElecs = data.badElecs;

% The same cleaning options apply to Granger, coherence, PPC, etc.
if add_params.linenoise_flag
    for iTrial = 1:numel(data.trial)
        data.trial{iTrial} = ft_preproc_dftfilter(data.trial{iTrial}, ...
            data.fsample, add_params.line_noise_freq, ...
            'dftreplace', add_params.ln_remove_method);
    end
end
if add_params.remove_bad_eletrod_flag
    goodChannels = true(numel(allLabels), 1);
    goodChannels(badElecs) = false;
    cfg = [];
    cfg.channel = allLabels(goodChannels);
    data = ft_selectdata(cfg, data);
end

% Compute the selected FC measure using the remaining channels.
[temp, fieldtripMethod] = getConnThisCondition(data, method, add_params);
freqVals = temp.freq;
% Keep the existing absolute-value convention, including for PPC.
fcValues = abs(temp.([fieldtripMethod 'spctrm']));
conn = nan(numel(allLabels), numel(allLabels), numel(freqVals));

if isfield(temp, 'labelcmb')
    % Bivariate GC returns pairs rather than a square channel matrix.
    pairLabels = regexprep(temp.labelcmb, '\[.*$', '');
    [~, sourceIdx] = ismember(pairLabels(:, 1), allLabels);
    [~, targetIdx] = ismember(pairLabels(:, 2), allLabels);
    if strcmp(fieldtripMethod, 'granger')
        % Bivariate self-GC is zero, as in the existing calculation.
        for iChannel = 1:numel(allLabels)
            conn(iChannel, iChannel, :) = 0;
        end
    end
    for iPair = 1:size(pairLabels, 1)
        conn(sourceIdx(iPair), targetIdx(iPair), :) = ...
            reshape(fcValues(iPair, :), 1, 1, []);
    end
else
    % Restore a square FC matrix to the original full channel order.
    [~, channelIdx] = ismember(temp.label, allLabels);
    conn(channelIdx, channelIdx, :) = fcValues;
end

% Mask bad output channels even when they were included in fitting, as before.
conn(badElecs, :, :) = NaN;
conn(:, badElecs, :) = NaN;
end

function [conn_stat, fieldtripMethod] = getConnThisCondition(data, method, add_params)

switch method
    case 'granger'
        switch add_params.granger_method
            case 'parametric'
                % First fit the AR model in time, then evaluate its spectra.
                cfgAR = [];
                cfgAR.method = add_params.mvar_method;
                cfgAR.feedback = 'none';
                cfgAR.order = add_params.mvar_order;
                if strcmp(cfgAR.method, 'biosig')
                    cfgAR.mvarmethod = add_params.mvar_algorithm;
                end
                cfgAR.keeptrials = 'no'; % Pooled fit; 'yes' is the placeholder above.
                cfgAR.demean = add_params.mvar_demean;
                cfgAR.zscore = add_params.mvar_zscore;
                cfgAR.taper = add_params.mvar_taper;

                cfgConn = [];
                cfgConn.method = 'granger';
                cfgConn.granger.conditional = 'no';
                switch add_params.granger_model
                    case 'bivariate'
                        % Fit each pair separately; this is the existing parametric path.
                        cfgAR.channelcmb = {'all', 'all'};
                        cfgConn.granger.sfmethod = 'bivariate';
                        cfgConn.channelcmb = {'all', 'all'};
                    case 'multivariate'
                        % Without channelcmb, FieldTrip fits all selected channels jointly.
                        cfgConn.granger.sfmethod = 'multivariate';
                end
                mdata = ft_mvaranalysis(cfgAR, data);
                cfgFreq = [];
                cfgFreq.foi = 0:add_params.mvar_freqstep:add_params.foilim_max;
                % The direct evaluator supports both pairwise and joint AR output.
                freqData = ft_freqanalysis_mvar(cfgFreq, mdata);
                conn_stat = ft_connectivityanalysis(cfgConn, freqData);

            case 'nonparametric'
                % Estimate Fourier coefficients, then factorize the pooled spectral matrix.
                freqData = getFourierSpectrum(data, add_params, 'fourier');
                cfgConn = [];
                cfgConn.method = 'granger';
                cfgConn.granger.conditional = 'no';
                switch add_params.granger_model
                    case 'bivariate'
                        % Wilson factorization is performed separately for each pair.
                        cfgConn.granger.sfmethod = 'bivariate';
                        cfgConn.channelcmb = {'all', 'all'};
                    case 'multivariate'
                        % Wilson factorization uses the joint spectral matrix.
                        cfgConn.granger.sfmethod = 'multivariate';
                end
                conn_stat = ft_connectivityanalysis(cfgConn, freqData);
        end
        fieldtripMethod = 'granger';

    case 'coh'
        % Average power and cross spectra, then compute coherence.
        freqData = getFourierSpectrum(data, add_params, 'powandcsd');
        freqData = ft_checkdata(freqData, 'cmbrepresentation', 'full');
        cfgConn = [];
        cfgConn.method = 'coh';
        conn_stat = ft_connectivityanalysis(cfgConn, freqData);
        fieldtripMethod = 'coh';

    case 'ppc'
        % PPC needs the individual Fourier observations, even for pooled FC.
        freqData = getFourierSpectrum(data, add_params, 'fourier');
        cfgConn = [];
        cfgConn.method = 'ppc';
        conn_stat = ft_connectivityanalysis(cfgConn, freqData);
        fieldtripMethod = 'ppc';

    otherwise
        % Keep the existing route for other Fourier-based FieldTrip methods.
        freqData = getFourierSpectrum(data, add_params, 'fourier');
        cfgConn = [];
        cfgConn.method = method;
        conn_stat = ft_connectivityanalysis(cfgConn, freqData);
        fieldtripMethod = method;
end
end

function freqData = getFourierSpectrum(data, add_params, outputType)

% The same Fourier settings are used by nonparametric Granger and phase FC.
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = add_params.taper;
if strcmp(cfg.taper, 'dpss')
    cfg.tapsmofrq = add_params.tapsmofrq;
end
cfg.output = outputType;
switch outputType
    case 'fourier'
        % Retain Fourier coefficients for pooling; this does not compute trial-wise GC.
        cfg.keeptrials = 'yes';
    case 'powandcsd'
        cfg.keeptrials = 'no';
end
cfg.foilim = [0 add_params.foilim_max];   % Keep the existing analyzed frequency range.
freqData = ft_freqanalysis(cfg, data);
end
