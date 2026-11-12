% Compare Granger spectra for different AR orders p.
% Detect both fs...Hz_p..._fmax...Hz and legacy p..._...Hz result folders.
% One figure per result root / sampling rate / line-noise setting; tabs are protocols.
% Average targets in each distance bin, then seeds, then subjects (omit NaNs).

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/home/venkatesh/Documents/IISc/data/Bsmart_Difffs_DiffOrderP
% /home/venkatesh/Documents/IISc/data/savedata_granger_fs_order_2026-11-05_18-08-03
%/home/venkatesh/Documents/IISc/data/Bsmart_Difffs_DiffOrderP
connDataRootList = { ...
    '/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_16-35-41-865/test1_fixed_order/biosig2'; ...
    '/home/venkatesh/Documents/IISc/data/OrderP_test'};
% Keep one entry above to plot just that result root.
mvarOrderList = [];                     % [] = all discovered orders; or [5 10 20 40 80].
subjects = {'053DR'};                   % One subject or e.g. {'053DR','019CKa'}.
% subjects = getGoodSubjectsBK1;         % Use all good subjects instead.
protocolNameList = {'G1'};%{'G1', 'G2', 'M1', 'M2'};
epochChoice = 'post';                   % 'pre', 'post', or 'combined' (mean of Pre/Post).
seeds = [17 18 48 16];                  % Oz, O2, POz, O1 in the saved 64-channel order.

% Frequency limits and sampling rates are read from the result folder names.
lineNoiseLabelList = {};                % {} = all; e.g. {'withoutLineNoise'}.
% Other labels: 'withLineNoise', 'noLineFrequenciesInRange'.
keeptrials = 'no';                      % This plot uses the pooled FC matrices.
FC_metric = 'granger_parametric';       % Or 'granger_parametric_multivariate'.
badEyeCondition = 'ep';
badTrialVersion = 'v8';

% Axis and line settings apply to every tab and panel.
xLimits = [];                          % [] = automatic; e.g. [0 100].
yLimits = [];                          % [] = automatic; e.g. [0 0.5].
gridSetting = 'on';
lineWidth = 1.5;

%%%%%%%%%%%%%%%%%%%%%%% Prepare distance groups %%%%%%%%%%%%%%%%%%%%%%%%%%
[electrodeGroups, binNames] = getElectrodeGroupsConn('rel', seeds, 'actiCap64_UOL');
numBins = numel(binNames);
for iSeed = 1:numel(seeds)
    for iBin = 1:numBins
        % Remove only this seed's self-connection from its target group.
        electrodeGroups{iSeed,iBin} = setdiff(electrodeGroups{iSeed,iBin}, seeds(iSeed), 'stable');
    end
end

% Discover available configurations instead of constructing one fixed folder name.
testsByRoot = cell(size(connDataRootList));
availableOrders = [];
for iRoot = 1:numel(connDataRootList)
    testsByRoot{iRoot} = findOrderTests(connDataRootList{iRoot}, keeptrials, lineNoiseLabelList);
    availableOrders = [availableOrders [testsByRoot{iRoot}.order]]; %#ok<AGROW>
end
if isempty(mvarOrderList)
    mvarOrderList = unique(availableOrders);
end
orderColors = lines(numel(mvarOrderList)); % Same p has the same color across all figures.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot each test group %%%%%%%%%%%%%%%%%%%%%%%%%
for iRoot = 1:numel(connDataRootList)
    connDataRoot = connDataRootList{iRoot};
    [~, rootName] = fileparts(connDataRoot);
    tests = testsByRoot{iRoot};
    tests = tests(ismember([tests.order], mvarOrderList));
    groupKeys = unique({tests.groupKey}, 'stable');
    for iGroup = 1:numel(groupKeys)
        groupTests = tests(strcmp({tests.groupKey}, groupKeys{iGroup}));
        [~, orderIdx] = sort([groupTests.order]);
        groupTests = groupTests(orderIdx);
        samplingLabel = getSamplingLabel(groupTests, subjects, protocolNameList, ...
            badEyeCondition, badTrialVersion, FC_metric);
        lineNoiseLabel = groupTests(1).noise;
        maxFrequency = groupTests(1).maxFrequency;
        fig = figure('Color', 'w', ...
            'Name', sprintf('%s | %s | fmax=%g Hz | %s', ...
                rootName, samplingLabel, maxFrequency, lineNoiseLabel), ...
            'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
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
                p = groupTests(iTest).order;
                iOrder = find(mvarOrderList == p, 1);
                subjectFC = [];
                plotFreq = [];
                numContributingSubjects = 0;
                for iSubject = 1:numel(subjects)
                    fileName = [protocolName '_' badEyeCondition '_' badTrialVersion '_' FC_metric '.mat'];
                    filePath = fullfile(groupTests(iTest).folder, subjects{iSubject}, fileName);
                    if ~isfile(filePath)
                        fprintf('Not saved yet / missing: %s\n', filePath);
                        continue
                    end
                    S = load(filePath, 'connPre', 'connPost', 'freqPre', 'freqPost', 'numGoodTrials');
                    if S.numGoodTrials == 0
                        continue
                    end
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
                    if isempty(C) || isempty(freq)
                        continue
                    end
                    if isempty(subjectFC)
                        plotFreq = freq;
                        subjectFC = nan(numel(subjects), numel(seeds), numBins, numel(freq));
                    elseif ~isequal(freq(:), plotFreq(:))
                        fprintf('Skipping mismatched subject frequency grid: %s\n', filePath);
                        continue
                    end

                    % Saved rows are sources; columns are targets. This is outgoing GC.
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
                if isempty(subjectFC) || numContributingSubjects == 0
                    continue
                end
                % First average seeds within each subject, then weight subjects equally.
                seedMean = mean(subjectFC, 2, 'omitnan');
                meanFC = reshape(mean(seedMean, 1, 'omitnan'), numBins, []);
                for iBin = 1:numBins
                    plot(ax(iBin), plotFreq, meanFC(iBin,:), 'Color', orderColors(iOrder,:), ...
                        'LineWidth', lineWidth, ...
                        'DisplayName', sprintf('p = %d (n = %d)', p, numContributingSubjects));
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
            title(tl, {sprintf('%s | %s | %s | %s', rootName, protocolName, epochChoice, samplingLabel), ...
                sprintf('fmax=%g Hz | %s | n = contributing subjects per curve', maxFrequency, lineNoiseLabel)}, ...
                'Interpreter', 'none');
        end
    end
end

function tests = findOrderTests(connDataRoot, keeptrials, lineNoiseLabelList)
tests = struct('folder', {}, 'order', {}, 'sampleRate', {}, 'maxFrequency', {}, 'noise', {}, 'groupKey', {});
folders = dir(connDataRoot);
for iFolder = 1:numel(folders)
    if ~folders(iFolder).isdir
        continue
    end
    name = folders(iFolder).name;
    token = regexp(name, ['^fs(?<fs>[0-9.eE+-]+)Hz_p(?<p>\d+)_fmax(?<fmax>[0-9.eE+-]+)Hz_' ...
        '(?:requested[0-9.eE+-]+Hz_)?(?<noise>withoutLineNoise|withLineNoise|noLineFrequenciesInRange)_keeptrials_(?<keep>.+)$'], 'names', 'once');
    if isempty(token)
        token = regexp(name, ['^p(?<p>\d+)_(?<fmax>[0-9.eE+-]+)Hz_' ...
            '(?<noise>withoutLineNoise|withLineNoise|noLineFrequenciesInRange)_keeptrials_(?<keep>.+)$'], 'names', 'once');
        if isempty(token)
            continue
        end
        sampleRate = NaN;               % Legacy names encode fmax, not sampling rate.
    else
        sampleRate = str2double(token.fs);
    end
    if ~strcmp(token.keep, keeptrials) || ...
            (~isempty(lineNoiseLabelList) && ~any(strcmp(token.noise, lineNoiseLabelList)))
        continue
    end
    test = struct;
    test.folder = fullfile(connDataRoot, name);
    test.order = str2double(token.p);
    test.sampleRate = sampleRate;
    test.maxFrequency = str2double(token.fmax);
    test.noise = token.noise;
    % Keep every non-order part of the folder name in the grouping key.
    test.groupKey = regexprep(name, '(^|_)p\d+_', '$1pORDER_');
    tests(end+1) = test; %#ok<AGROW>
end
end

function samplingLabel = getSamplingLabel(tests, subjects, protocolNameList, badEyeCondition, badTrialVersion, FC_metric)
if ~isnan(tests(1).sampleRate)
    samplingLabel = sprintf('fs=%g Hz', tests(1).sampleRate);
    return
end
samplingLabel = 'sampling rate not recorded';
for iTest = 1:numel(tests)
    for iSubject = 1:numel(subjects)
        for iProtocol = 1:numel(protocolNameList)
            fileName = [protocolNameList{iProtocol} '_' badEyeCondition '_' badTrialVersion '_' FC_metric '.mat'];
            filePath = fullfile(tests(iTest).folder, subjects{iSubject}, fileName);
            if ~isfile(filePath), continue; end
            meta = load(filePath, 'add_params');
            if isfield(meta, 'add_params') && isfield(meta.add_params, 'resamplefs')
                if isempty(meta.add_params.resamplefs)
                    samplingLabel = 'native sampling (no resampling)';
                else
                    samplingLabel = sprintf('fs=%g Hz', meta.add_params.resamplefs);
                end
                return
            end
        end
    end
end
end
