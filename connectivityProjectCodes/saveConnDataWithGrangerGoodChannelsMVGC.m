function saveConnDataWithGrangerGoodChannelsMVGC(subjectName, protocolNameList, badEyeCondition, badTrialVersion, ftDataFolder, connMethod, stRange, connDataFolder, add_params)
% MVGC1 adapter for the existing FieldTrip input files and Pre/Post windows.
% granger_model='bivariate': fit each pair, nonconditional GC.
% granger_model='multivariate', mvgc_conditional='yes': joint, conditional GC.
% granger_model='multivariate', mvgc_conditional='no': joint fit, marginal pair GC.
% Joint marginal GC uses autocovariance submatrices, regardless of spectral_method.
% Saved matrices use SOURCE x TARGET x frequency; MVGC's native axes are reversed.

initialiseMVGC(add_params.mvgc_root);
if ~strcmp(add_params.keeptrials, 'no')
    error('MVGC adapter supports pooled trials only: keeptrials = ''no''.');
end
for iProtocol = 1:numel(protocolNameList)
    protocolName = protocolNameList{iProtocol};
    saveProtocol(subjectName, protocolName, badEyeCondition, badTrialVersion, ...
        ftDataFolder, connMethod, stRange, connDataFolder, add_params);
end
end

function saveProtocol(subjectName, protocolName, badEyeCondition, badTrialVersion, ftDataFolder, connMethod, stRange, connDataFolder, add_params)
subjectFolder = fullfile(connDataFolder, subjectName);
if ~isfolder(subjectFolder), mkdir(subjectFolder); end
ftDataFileName = fullfile(ftDataFolder, subjectName, ...
    [protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
S = load(ftDataFileName, 'data', 'numGoodTrials');
numGoodTrials = S.numGoodTrials;
fsUsed = [];
connPre = []; connPost = []; freqPre = []; freqPost = [];
gcTimePre = []; gcTimePost = [];
diagnosticsPre = struct; diagnosticsPost = struct;
if numGoodTrials > 0
    data = S.data;
    if ~isempty(add_params.resamplefs)
        if add_params.resamplefs > data.fsample
            error('Requested sampling rate exceeds input rate in %s.', ftDataFileName);
        elseif add_params.resamplefs < data.fsample
            badElecs = data.badElecs;
            cfg = []; cfg.resamplefs = add_params.resamplefs;
            cfg.method = 'resample'; cfg.detrend = 'no';
            data = ft_resampledata(cfg, data);
            data.badElecs = badElecs;
        end
    end
    fsUsed = data.fsample;
    add_params.line_noise_freq = add_params.line_noise_freq(add_params.line_noise_freq < fsUsed/2);
    add_params.linenoise_flag = add_params.linenoise_flag && ~isempty(add_params.line_noise_freq);
    cfg = []; cfg.toilim = [-diff(stRange) + 1/fsUsed, 0];
    dataPre = ft_redefinetrial(cfg, data);
    cfg.toilim = [stRange(1) + 1/fsUsed, stRange(2)];
    dataPost = ft_redefinetrial(cfg, data);
    [connPre, freqPre, gcTimePre, diagnosticsPre] = calculateGC(dataPre, add_params);
    [connPost, freqPost, gcTimePost, diagnosticsPost] = calculateGC(dataPost, add_params);
    if diagnosticsPre.numFailedModels > 0 || diagnosticsPost.numFailedModels > 0
        fprintf('MVGC %s %s: failed models Pre=%d, Post=%d; inspect saved diagnostics.\n', ...
            subjectName, protocolName, diagnosticsPre.numFailedModels, diagnosticsPost.numFailedModels);
    end
end
connDimord = 'source_target_freq';
gcTimeDimord = 'source_target';
mvgcSourcePath = fileparts(fileparts(which('tsdata_to_var')));
global mvgc_version
mvgcToolboxVersion = mvgc_version;
outputFile = fullfile(subjectFolder, ...
    [protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']);
save(outputFile, 'connPre', 'connPost', 'freqPre', 'freqPost', ...
    'gcTimePre', 'gcTimePost', 'diagnosticsPre', 'diagnosticsPost', ...
    'numGoodTrials', 'ftDataFileName', 'add_params', 'fsUsed', 'connDimord', ...
    'gcTimeDimord', 'mvgcSourcePath', 'mvgcToolboxVersion', '-v7.3');
end

function [conn, freqVals, gcTime, diagnostics] = calculateGC(data, options)
numAllChannels = numel(data.label);
fitIndices = 1:numAllChannels;
if options.remove_bad_eletrod_flag
    fitIndices = setdiff(fitIndices, data.badElecs, 'stable');
end
numChannels = numel(fitIndices);
lengths = cellfun(@(x) size(x,2), data.trial);
if isempty(lengths) || any(lengths ~= lengths(1))
    error('MVGC needs nonempty, equal-length trials within each Pre/Post window.');
end
X = zeros(numChannels, lengths(1), numel(data.trial));
for iTrial = 1:numel(data.trial)
    trial = double(data.trial{iTrial}(fitIndices,:));
    if options.linenoise_flag
        trial = ft_preproc_dftfilter(trial, data.fsample, options.line_noise_freq, ...
            'dftreplace', options.ln_remove_method);
    end
    if strcmp(options.mvgc_demean, 'pertrial')
        trial = trial - mean(trial, 2);
    elseif ~strcmp(options.mvgc_demean, 'ensemble')
        error('mvgc_demean must be ''ensemble'' or ''pertrial''.');
    end
    X(:,:,iTrial) = trial;
end
% MVGC always removes the pooled channel mean, without joining trial boundaries.
if strcmp(options.mvar_zscore, 'yes')
    scale = std(reshape(X,numChannels,[]), 0, 2);
    X = X ./ reshape(scale,[],1,1);
end
if isempty(options.mvgc_fres)
    fres = max(1, ceil((data.fsample/2)/options.mvar_freqstep));
else
    fres = options.mvgc_fres;
end
fullFreq = sfreqs(fres, data.fsample)';
keepFreq = fullFreq <= min(options.foilim_max, data.fsample/2) + 10*eps(data.fsample);
freqVals = fullFreq(keepFreq);
fullGC = nan(numChannels,numChannels,fres+1);
fullTime = nan(numChannels);
diagnostics = struct('fittedChannelIndices', fitIndices, ...
    'fittedChannelLabels', {data.label(fitIndices)}, 'sampleRate', data.fsample, ...
    'samplesPerTrial', lengths(1), 'numTrials', numel(data.trial), ...
    'frequencyResolution', fres, 'actualFrequencyStepHz', data.fsample/(2*fres), ...
    'centering', options.mvgc_demean, 'models', {{}});
if strcmp(options.granger_model, 'bivariate')
    if strcmp(options.mvgc_conditional, 'yes')
        error('Conditional GC needs granger_model = ''multivariate''.');
    end
    for i = 1:numChannels-1
        for j = i+1:numChannels
            pair = [i j];
            [f,F,info] = fitAndEvaluate(X(pair,:,:), options, fres);
            info.channelIndices = fitIndices(pair);
            diagnostics.models{end+1} = info;
            fullGC(pair,pair,:) = f;
            fullTime(pair,pair) = F;
        end
    end
elseif strcmp(options.granger_model, 'multivariate')
    [fullGC,fullTime,info] = fitAndEvaluate(X, options, fres);
    info.channelIndices = fitIndices;
    diagnostics.models{1} = info;
else
    error('granger_model must be ''bivariate'' or ''multivariate''.');
end
for i = 1:numChannels
    fullGC(i,i,:) = 0; fullTime(i,i) = 0;
end
conn = nan(numAllChannels,numAllChannels,nnz(keepFreq));
gcTime = nan(numAllChannels);
% MVGC f(target,source,freq) -> existing saved format conn(source,target,freq).
conn(fitIndices,fitIndices,:) = permute(fullGC(:,:,keepFreq),[2 1 3]);
gcTime(fitIndices,fitIndices) = fullTime';
conn(data.badElecs,:,:) = NaN; conn(:,data.badElecs,:) = NaN;
gcTime(data.badElecs,:) = NaN; gcTime(:,data.badElecs) = NaN;
diagnostics.numFailedModels = sum(cellfun(@(d) ~d.success, diagnostics.models));
end

function [f,F,info] = fitAndEvaluate(X, options, fres)
n = size(X,1);
f = nan(n,n,fres+1); F = nan(n);
info = struct('order', options.mvar_order, 'regmode', options.mvgc_regmode, ...
    'spectralMethod', options.mvgc_spectral_method, 'success', false, ...
    'message', '', 'spectralRadius', NaN, 'autocovarianceInfo', struct, ...
    'spectralIntegralMaxAbsError', NaN);
if options.mvar_order >= size(X,2)
    info.message = 'Model order must be smaller than samples per trial.'; return
end
try
    [A,V] = tsdata_to_var(X, options.mvar_order, options.mvgc_regmode);
    if any(~isfinite(A(:))) || any(~isfinite(V(:)))
        info.message = 'VAR fit returned nonfinite coefficients/covariance.'; return
    end
    info.spectralRadius = var_specrad(A);
    [~,cholFlag] = chol(V);
    if info.spectralRadius >= 1 || cholFlag ~= 0
        info.message = 'Unstable VAR or non-positive-definite residual covariance.'; return
    end
    jointMarginal = n > 2 && strcmp(options.mvgc_conditional, 'no');
    if jointMarginal || strcmp(options.mvgc_spectral_method, 'autocov')
        % MVGC1's Schur path calls an internal MATLAB function removed in
        % newer releases. Use MVGC's own iterative solver when unavailable.
        useIterative = exist('lyapslv','file') == 0;
        [G,acInfo] = var_to_autocov(A,V,options.mvgc_acmaxlags,options.mvgc_acdectol,useIterative);
        if useIterative
            acInfo.lyapunovSolver = 'accelerated iterative';
        else
            acInfo.lyapunovSolver = 'Schur';
        end
        info.autocovarianceInfo = acInfo;
        info.spectralMethod = 'autocov';
        if acInfo.error ~= 0
            info.message = acInfo.errmsg; return
        end
        if jointMarginal
            for i = 1:n-1
                for j = i+1:n
                    pair = [i j];
                    f(pair,pair,:) = autocov_to_spwcgc(G(pair,pair,:),fres);
                    F(pair,pair) = autocov_to_pwcgc(G(pair,pair,:));
                end
            end
        else
            f = autocov_to_spwcgc(G,fres);
            F = autocov_to_pwcgc(G);
        end
    elseif strcmp(options.mvgc_spectral_method, 'statespace')
        f = var_to_spwcgc(A,V,fres);
        F = var_to_pwcgc(A,V);
    else
        error('mvgc_spectral_method must be ''statespace'' or ''autocov''.');
    end
    offDiagonal = ~eye(n);
    flatF = reshape(f,n*n,[]);
    if ~isreal(f) || any(~isfinite(flatF(offDiagonal(:),:)), 'all') ...
            || any(~isfinite(F(offDiagonal)))
        info.message = 'GC calculation failed; inspect the model and MVGC warnings.';
        f(:) = NaN; F(:) = NaN; return
    end
    % Check the complete Nyquist spectrum, before trimming to foilim_max.
    integrated = smvgc_to_mvgc(f);
    info.spectralIntegralMaxAbsError = max(abs(integrated(offDiagonal)-F(offDiagonal)));
    info.success = true;
catch modelError
    info.message = modelError.message;
    f(:) = NaN; F(:) = NaN;
end
% No abs() or clipping: invalid/negative GC must not be disguised as valid GC.
end

function initialiseMVGC(mvgcRoot)
if exist('tsdata_to_var','file') == 2 && exist('var_to_spwcgc','file') == 2
    return
end
warningState = warning; rngState = rng;
restore = onCleanup(@() restoreSession(warningState,rngState));
run(fullfile(mvgcRoot,'startup.m')); % Official startup selects paths; never use genpath.
end

function restoreSession(warningState,rngState)
warning(warningState); rng(rngState);
end
