% Compare saved Granger spectra across sampling rates at each fixed AR order.
% Like Plot_versusOrderP: protocol tabs, distance-bin panels, outgoing GC.
% Average targets in each bin, then seeds, then subjects (omit NaNs).
% This script only reads saved results; it does not recompute connectivity.
% Optional PPC overlay uses a dashed black curve and a separate right y-axis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_18-57-55-545/test1_fixed_order/biosig2
%- p60
%/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_16-35-41-865/test1_fixed_order/biosig2
%- p20
connDataRoot = '/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_18-57-55-545/test1_fixed_order/biosig2';
sampleRateList = [];                    % [] = all saved rates; e.g. [50 100 200 500].
mvarOrderList = [];                     % [] = all saved orders, one figure per p.
subjects = {'053DR'};
% subjects = getGoodSubjectsBK1;         % Use all good subjects instead.
protocolNameList = {'G2', 'M2'};               % Or {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'}.
epochChoice = 'post';                   % 'pre', 'post', or 'combined' (mean of Pre/Post).
seeds = [17 18 48 16];                  % Oz, O2, POz, O1 in the saved channel order.

overlayPPC = true;
ppcDataRoot = '/home/venkatesh/Documents/IISc/data/PPC_LineNoiseRemoved';
ppcYLimits = [];                       % [] = automatic right axis; e.g. [0 1].

keeptrials = 'no';
FC_metric = 'granger_parametric';
badEyeCondition = 'ep';
badTrialVersion = 'v8';
xLimits = [];                          % [] = full saved frequency range; or [0 100].
yLimits = [];                          % [] = automatic; e.g. [0 0.5].
gridSetting = 'on';
lineWidth = 1.5;

%%%%%%%%%%%%%%%%%%%%%%% Prepare distance groups %%%%%%%%%%%%%%%%%%%%%%%%%%
[electrodeGroups, binNames] = getElectrodeGroupsConn('rel', seeds, 'actiCap64_UOL');
numBins = numel(binNames);
for iSeed = 1:numel(seeds)
    for iBin = 1:numBins
        electrodeGroups{iSeed,iBin} = setdiff(electrodeGroups{iSeed,iBin}, seeds(iSeed), 'stable');
    end
end

tests = findSampleRateTests(connDataRoot, keeptrials);
if isempty(tests)
    fprintf('No sampling-rate result folders found in %s\n', connDataRoot);
    return
end
if isempty(sampleRateList), sampleRateList = unique([tests.sampleRate]); end
if isempty(mvarOrderList), mvarOrderList = unique([tests.order]); end
sampleRateList = sort(unique(sampleRateList));
tests = tests(ismember([tests.sampleRate], sampleRateList) & ismember([tests.order], mvarOrderList));
rateColors = turbo(numel(sampleRateList)); % Distinct colors even with more than seven rates.
groupKeys = unique({tests.groupKey}, 'stable');
[~, rootName] = fileparts(connDataRoot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot each fixed p %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iGroup = 1:numel(groupKeys)
    groupTests = tests(strcmp({tests.groupKey}, groupKeys{iGroup}));
    [~, rateIdx] = sort([groupTests.sampleRate]);
    groupTests = groupTests(rateIdx);
    p = groupTests(1).order;
    noiseLabel = groupTests(1).noisePolicy;
    fig = figure('Color', 'w', 'NumberTitle', 'off', ...
        'Name', sprintf('%s | p=%d | %s', rootName, p, noiseLabel), ...
        'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tabGroup = uitabgroup(fig);

    for iProtocol = 1:numel(protocolNameList)
        protocolName = protocolNameList{iProtocol};
        tab = uitab(tabGroup, 'Title', protocolName);
        tl = tiledlayout(tab, ceil(numBins/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax = gobjects(1, numBins);
        for iBin = 1:numBins
            ax(iBin) = nexttile(tl);
            hold(ax(iBin), 'on');
        end

        for iTest = 1:numel(groupTests)
            test = groupTests(iTest);
            iRate = find(sampleRateList == test.sampleRate, 1);
            subjectFC = [];
            plotFreq = [];
            numContributingSubjects = 0;
            for iSubject = 1:numel(subjects)
                fileName = [protocolName '_' badEyeCondition '_' badTrialVersion '_' FC_metric '.mat'];
                filePath = fullfile(test.folder, subjects{iSubject}, fileName);
                if ~isfile(filePath)
                    fprintf('Not saved yet / missing: %s\n', filePath);
                    continue
                end
                S = load(filePath, 'connPre', 'connPost', 'freqPre', 'freqPost', 'numGoodTrials');
                if S.numGoodTrials == 0, continue; end
                switch epochChoice
                    case 'pre'
                        C = S.connPre; freq = S.freqPre;
                    case 'post'
                        C = S.connPost; freq = S.freqPost;
                    case 'combined'
                        if ~isequal(S.freqPre(:), S.freqPost(:))
                            fprintf('Skipping mismatched Pre/Post frequency grids: %s\n', filePath);
                            continue
                        end
                        C = (S.connPre + S.connPost)/2; freq = S.freqPost;
                end
                if isempty(C) || isempty(freq), continue; end
                if isempty(subjectFC)
                    plotFreq = freq;
                    subjectFC = nan(numel(subjects), numel(seeds), numBins, numel(freq));
                elseif ~isequal(freq(:), plotFreq(:))
                    fprintf('Skipping mismatched subject frequency grid: %s\n', filePath);
                    continue
                end
                % Saved rows are sources and columns are targets: outgoing GC.
                for iSeed = 1:numel(seeds)
                    for iBin = 1:numBins
                        targets = electrodeGroups{iSeed,iBin};
                        subjectFC(iSubject,iSeed,iBin,:) = reshape( ...
                            mean(C(seeds(iSeed),targets,:), 2, 'omitnan'), 1, 1, 1, []);
                    end
                end
                if any(isfinite(subjectFC(iSubject,:,:,:)), 'all')
                    numContributingSubjects = numContributingSubjects + 1;
                end
            end
            if isempty(subjectFC) || numContributingSubjects == 0, continue; end
            seedMean = mean(subjectFC, 2, 'omitnan');
            meanFC = reshape(mean(seedMean, 1, 'omitnan'), numBins, []);
            for iBin = 1:numBins
                % Use each rate's actual grid; do not interpolate or extrapolate.
                plot(ax(iBin), plotFreq, meanFC(iBin,:), ...
                    'Color', rateColors(iRate,:), 'LineWidth', lineWidth, ...
                    'DisplayName', sprintf('Fs = %g Hz (n = %d)', test.sampleRate, numContributingSubjects), ...
                    'UserData', struct('sampleRate',test.sampleRate,'order',p,'folder',test.folder));
            end
            fprintf('Plotted %s | p=%d | Fs=%g Hz | fmax=%g Hz | %d subject(s)\n', ...
                protocolName, p, test.sampleRate, max(plotFreq), numContributingSubjects);
        end

        if overlayPPC
            [ppcFreq, meanPPC, numPPCSubjects] = loadPPCMean(ppcDataRoot, subjects, ...
                protocolName, badEyeCondition, badTrialVersion, epochChoice, seeds, electrodeGroups);
            if numPPCSubjects > 0
                for iBin = 1:numBins
                    yyaxis(ax(iBin), 'right');
                    plot(ax(iBin), ppcFreq, meanPPC(iBin,:), 'k--', 'LineWidth', 2, ...
                        'DisplayName', sprintf('PPC (n = %d, right axis)', numPPCSubjects), ...
                        'Tag', 'PPCOverlay', 'UserData', struct('metric','ppc','folder',ppcDataRoot));
                    ylabel(ax(iBin), 'Seed-averaged PPC');
                    ax(iBin).YAxis(2).Color = [0 0 0];
                    if ~isempty(ppcYLimits), ylim(ax(iBin), ppcYLimits); end
                    yyaxis(ax(iBin), 'left');
                    ax(iBin).YAxis(1).Color = [0 0 0];
                end
                fprintf('Overlaid PPC | %s | %s | %d subject(s) | saved grid %g-%g Hz\n', ...
                    protocolName, epochChoice, numPPCSubjects, min(ppcFreq), max(ppcFreq));
            end
        end

        for iBin = 1:numBins
            title(ax(iBin), binNames{iBin});
            xlabel(ax(iBin), 'Frequency (Hz)');
            ylabel(ax(iBin), 'Seed-averaged outgoing GC');
            grid(ax(iBin), gridSetting);
            if ~isempty(findobj(ax(iBin), 'Type', 'line'))
                legend(ax(iBin), 'show', 'Location', 'best');
            end
            if ~isempty(xLimits), xlim(ax(iBin), xLimits); end
            if ~isempty(yLimits), ylim(ax(iBin), yLimits); end
        end
        title(tl, {sprintf('%s | %s | %s | fixed p=%d', rootName, protocolName, epochChoice, p), ...
            sprintf('%s | each Fs uses its saved frequency range | n = contributing subjects', noiseLabel)}, ...
            'Interpreter', 'none');
    end
end

function [plotFreq, meanFC, numContributingSubjects] = loadPPCMean(root, subjects, protocol, badEyeCondition, badTrialVersion, epochChoice, seeds, electrodeGroups)
% Apply the same target -> seed -> subject averaging to the saved PPC values.
% PPC is nondirectional. Preserve its native values and native frequency grid.
numBins = size(electrodeGroups,2);
plotFreq = [];
meanFC = [];
subjectFC = [];
numContributingSubjects = 0;
for iSubject = 1:numel(subjects)
    fileName = [protocol '_' badEyeCondition '_' badTrialVersion '_ppc.mat'];
    filePath = fullfile(root, subjects{iSubject}, fileName);
    if ~isfile(filePath)
        fprintf('PPC missing: %s\n', filePath);
        continue
    end
    S = load(filePath, 'connPre', 'connPost', 'freqPre', 'freqPost', 'numGoodTrials');
    if S.numGoodTrials == 0, continue; end
    switch epochChoice
        case 'pre'
            C = S.connPre; freq = S.freqPre;
        case 'post'
            C = S.connPost; freq = S.freqPost;
        case 'combined'
            if ~isequal(S.freqPre(:), S.freqPost(:))
                fprintf('Skipping mismatched PPC Pre/Post frequency grids: %s\n', filePath);
                continue
            end
            C = (S.connPre + S.connPost)/2; freq = S.freqPost;
    end
    if isempty(C) || isempty(freq), continue; end
    if isempty(subjectFC)
        plotFreq = freq;
        subjectFC = nan(numel(subjects), numel(seeds), numBins, numel(freq));
    elseif ~isequal(freq(:), plotFreq(:))
        fprintf('Skipping mismatched PPC subject frequency grid: %s\n', filePath);
        continue
    end
    for iSeed = 1:numel(seeds)
        for iBin = 1:numBins
            targets = electrodeGroups{iSeed,iBin};
            subjectFC(iSubject,iSeed,iBin,:) = reshape( ...
                mean(C(seeds(iSeed),targets,:), 2, 'omitnan'), 1, 1, 1, []);
        end
    end
    if any(isfinite(subjectFC(iSubject,:,:,:)), 'all')
        numContributingSubjects = numContributingSubjects + 1;
    end
end
if numContributingSubjects > 0
    meanFC = reshape(mean(mean(subjectFC, 2, 'omitnan'), 1, 'omitnan'), numBins, []);
end
end

function tests = findSampleRateTests(root, keeptrials)
tests = struct('folder',{},'sampleRate',{},'order',{},'noisePolicy',{},'groupKey',{});
folders = dir(root);
for iFolder = 1:numel(folders)
    if ~folders(iFolder).isdir, continue; end
    token = regexp(folders(iFolder).name, ...
        ['^fs(?<fs>[0-9.eE+-]+)Hz_p(?<p>\d+)_fmax(?<fmax>[0-9.eE+-]+)Hz_' ...
        '(?:requested[0-9.eE+-]+Hz_)?(?<noise>withoutLineNoise|withLineNoise|noLineFrequenciesInRange)_keeptrials_(?<keep>.+)$'], ...
        'names', 'once');
    if isempty(token) || ~strcmp(token.keep, keeptrials), continue; end
    % At Fs=50/100 no configured 50/100-Hz line frequencies are in range.
    % Keep those curves with the other rates from the same removal policy.
    if strcmp(token.noise, 'withLineNoise')
        noisePolicy = 'line noise retained';
    else
        noisePolicy = 'line noise removed where in range';
    end
    test = struct('folder',fullfile(root,folders(iFolder).name), ...
        'sampleRate',str2double(token.fs),'order',str2double(token.p), ...
        'noisePolicy',noisePolicy,'groupKey',[token.p '_' noisePolicy]);
    tests(end+1) = test; %#ok<AGROW>
end
end
