% function figures = Plot_versusOrderP_FixedFs(sampleRate, subjects, protocolNameList)
% Overlay saved AR orders at one Fs, with optional PPC on the right y-axis.
% Examples (each call opens a new figure and keeps earlier figures):
%   Plot_versusOrderP_FixedFs(200, '053DR', 'G1')
%   Plot_versusOrderP_FixedFs(500, '053DR', 'G1')
% Subjects/protocols may also be cell arrays. Missing files are skipped.
% Average targets within distance bins, then seeds, then subjects (omit NaNs).

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_versusOrderP_FixedFs(200, '', 'G1');
% Plot_versusOrderP_FixedFs(500, '053DR', 'G1');
sampleRate = 500;
subjects = '053DR';
protocolNameList = 'G1';
% if nargin < 1 || isempty(sampleRate), sampleRate = 200; end
% if nargin < 2 || isempty(subjects), subjects = {'053DR'}; end
% if nargin < 3 || isempty(protocolNameList), protocolNameList = {'G1'}; end
subjects = cellstr(subjects);
protocolNameList = cellstr(protocolNameList);
%/home/venkatesh/Documents/IISc/data/savedata_granger_fs_order_2026-11-05_18-08-03
%/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_16-35-41-865/test1_fixed_order/biosig2
connDataRoot = '/home/venkatesh/Documents/IISc/data/savedata_biosig_parallel_sweeps_2026-11-06_16-35-41-865/test2_order_sweep/biosig2';
mvarOrderList = [];                     % [] = all discovered orders; e.g. [5 10 25 50].
epochChoice = 'post';                   % 'pre', 'post', or 'combined'.
seeds = [17 18 48 16];                  % Oz, O2, POz, O1.
overlayPPC = true;
ppcDataRoot = '/home/venkatesh/Documents/IISc/data/PPC_LineNoiseRemoved';
ppcYLimits = [];                        % [] = automatic right axis.
keeptrials = 'no';
FC_metric = 'granger_parametric';
badEyeCondition = 'ep';
badTrialVersion = 'v8';
xLimits = [];                          % [] = [0 Fs/2]; or e.g. [0 100].
yLimits = [];                          % [] = automatic GC axis.
gridSetting = 'on';
lineWidth = 1.5;

figures = gobjects(0);
tests = findOrderSweepTests(connDataRoot, keeptrials);
if ~isempty(mvarOrderList)
    tests = tests(ismember([tests.order], mvarOrderList));
end
allOrders = unique([tests.order]);    % Keep a given p the same color at 200/500 Hz.
tests = tests([tests.sampleRate] == sampleRate);
if isempty(tests)
    fprintf('No matching order-sweep folders for Fs=%g Hz in %s\n', sampleRate, connDataRoot);
    return
end
orderColors = turbo(numel(allOrders));
[electrodeGroups, binNames] = getElectrodeGroupsConn('rel', seeds, 'actiCap64_UOL');
numBins = numel(binNames);
for iSeed = 1:numel(seeds)
    for iBin = 1:numBins
        electrodeGroups{iSeed,iBin} = setdiff(electrodeGroups{iSeed,iBin}, seeds(iSeed), 'stable');
    end
end
groupKeys = unique({tests.groupKey}, 'stable');
[~, rootName] = fileparts(connDataRoot);

for iGroup = 1:numel(groupKeys)
    groupTests = tests(strcmp({tests.groupKey}, groupKeys{iGroup}));
    [~, orderIdx] = sort([groupTests.order]);
    groupTests = groupTests(orderIdx);
    fig = gobjects(0);
    for iProtocol = 1:numel(protocolNameList)
        protocol = protocolNameList{iProtocol};
        curves = struct('freq',{},'meanFC',{},'n',{},'order',{},'folder',{});
        for iTest = 1:numel(groupTests)
            test = groupTests(iTest);
            [freq, meanFC, n] = readMeanSpectrum(test.folder, subjects, protocol, FC_metric, ...
                badEyeCondition, badTrialVersion, epochChoice, seeds, electrodeGroups);
            if n == 0, continue; end
            curves(end+1) = struct('freq',freq,'meanFC',meanFC,'n',n, ...
                'order',test.order,'folder',test.folder); %#ok<AGROW>
        end
        ppcFreq = []; meanPPC = []; nPPC = 0;
        if overlayPPC
            [ppcFreq, meanPPC, nPPC] = readMeanSpectrum(ppcDataRoot, subjects, protocol, 'ppc', ...
                badEyeCondition, badTrialVersion, epochChoice, seeds, electrodeGroups);
        end
        if isempty(curves) && nPPC == 0
            fprintf('No usable GC or PPC results: Fs=%g Hz | %s. No tab created.\n', sampleRate, protocol);
            continue
        end
        % Create figures/tabs only when there is something to plot.
        if isempty(fig)
            fig = figure('Color','w','NumberTitle','off', ...
                'Name',sprintf('%s | Fs=%g Hz | orders + PPC | %s',rootName,sampleRate,groupTests(1).noisePolicy), ...
                'Units','normalized','OuterPosition',[0 0 1 1]);
            figures(end+1) = fig; %#ok<AGROW>
            tabGroup = uitabgroup(fig);
        end
        tab = uitab(tabGroup,'Title',protocol);
        tl = tiledlayout(tab,ceil(numBins/2),2,'TileSpacing','compact','Padding','compact');
        for iBin = 1:numBins
            ax = nexttile(tl);
            hold(ax,'on');
            for iCurve = 1:numel(curves)
                curve = curves(iCurve);
                iColor = find(allOrders == curve.order,1);
                plot(ax,curve.freq,curve.meanFC(iBin,:), ...
                    'Color',orderColors(iColor,:),'LineWidth',lineWidth, ...
                    'DisplayName',sprintf('p = %d (n = %d)',curve.order,curve.n), ...
                    'Tag','GCOrderCurve','UserData',struct('sampleRate',sampleRate,'order',curve.order,'folder',curve.folder));
            end
            ylabel(ax,'Seed-averaged outgoing GC');
            if ~isempty(yLimits), ylim(ax,yLimits); end
            if nPPC > 0
                yyaxis(ax,'right');
                plot(ax,ppcFreq,meanPPC(iBin,:),'k--','LineWidth',2, ...
                    'DisplayName',sprintf('PPC (n = %d, right axis)',nPPC), ...
                    'Tag','PPCOverlay','UserData',struct('metric','ppc','folder',ppcDataRoot));
                ylabel(ax,'Seed-averaged PPC');
                ax.YAxis(2).Color = [0 0 0];
                if ~isempty(ppcYLimits), ylim(ax,ppcYLimits); end
                yyaxis(ax,'left');
                ax.YAxis(1).Color = [0 0 0];
            end
            title(ax,binNames{iBin});
            xlabel(ax,'Frequency (Hz)');
            grid(ax,gridSetting);
            legend(ax,'show','Location','best');
            if isempty(xLimits), xlim(ax,[0 sampleRate/2]); else, xlim(ax,xLimits); end
        end
        title(tl,{sprintf('%s | Fs=%g Hz | %s | %s | %s', ...
            rootName,sampleRate,strjoin(subjects,', '),protocol,epochChoice), ...
            sprintf('%s | n = contributing subjects per curve',groupTests(1).noisePolicy)}, ...
            'Interpreter','none');
        fprintf('Plotted Fs=%g Hz | %s | p=%s | PPC subjects=%d\n', ...
            sampleRate,protocol,mat2str([curves.order]),nPPC);
    end
end
% end

function [plotFreq, meanFC, n] = readMeanSpectrum(root, subjects, protocol, metric, badEyeCondition, badTrialVersion, epochChoice, seeds, electrodeGroups)
plotFreq = []; meanFC = []; subjectFC = []; n = 0;
numBins = size(electrodeGroups,2);
required = {'connPre','connPost','freqPre','freqPost','numGoodTrials'};
for iSubject = 1:numel(subjects)
    filePath = fullfile(root,subjects{iSubject},[protocol '_' badEyeCondition '_' badTrialVersion '_' metric '.mat']);
    if ~isfile(filePath)
        fprintf('Missing / not saved yet, skipping: %s\n',filePath);
        continue
    end
    try
        S = load(filePath,required{:});
    catch readError
        fprintf('Cannot read yet, skipping: %s (%s)\n',filePath,readError.message);
        continue
    end
    if ~all(isfield(S,required))
        fprintf('Incomplete saved fields, skipping: %s\n',filePath);
        continue
    end
    if isempty(S.numGoodTrials) || S.numGoodTrials == 0, continue; end
    switch epochChoice
        case 'pre'
            C = S.connPre; freq = S.freqPre;
        case 'post'
            C = S.connPost; freq = S.freqPost;
        case 'combined'
            if ~isequal(S.freqPre(:),S.freqPost(:))
                fprintf('Mismatched Pre/Post frequency grids, skipping: %s\n',filePath);
                continue
            end
            C = (S.connPre+S.connPost)/2; freq = S.freqPost;
    end
    if isempty(C) || isempty(freq), continue; end
    if isempty(subjectFC)
        plotFreq = freq;
        subjectFC = nan(numel(subjects),numBins,numel(freq));
    elseif ~isequal(freq(:),plotFreq(:))
        fprintf('Mismatched subject frequency grids, skipping: %s\n',filePath);
        continue
    end
    seedFC = nan(numel(seeds),numBins,numel(freq));
    for iSeed = 1:numel(seeds)
        for iBin = 1:numBins
            targets = electrodeGroups{iSeed,iBin};
            seedFC(iSeed,iBin,:) = mean(C(seeds(iSeed),targets,:),2,'omitnan');
        end
    end
    subjectFC(iSubject,:,:) = mean(seedFC,1,'omitnan');
    if any(isfinite(subjectFC(iSubject,:,:)),'all'), n = n+1; end
end
if n > 0
    meanFC = reshape(mean(subjectFC,1,'omitnan'),numBins,[]);
end
end

function tests = findOrderSweepTests(root, keeptrials)
tests = struct('folder',{},'sampleRate',{},'order',{},'noisePolicy',{},'groupKey',{});
folders = dir(root);
for iFolder = 1:numel(folders)
    if ~folders(iFolder).isdir, continue; end
    token = regexp(folders(iFolder).name, ...
        ['^fs(?<fs>[0-9.eE+-]+)Hz_p(?<p>\d+)_fmax(?<fmax>[0-9.eE+-]+)Hz_' ...
        '(?:requested[0-9.eE+-]+Hz_)?(?<noise>withoutLineNoise|withLineNoise|noLineFrequenciesInRange)_keeptrials_(?<keep>.+)$'], ...
        'names','once');
    if isempty(token) || ~strcmp(token.keep,keeptrials), continue; end
    if strcmp(token.noise,'withLineNoise')
        noisePolicy = 'line noise retained';
    else
        noisePolicy = 'line noise removed where in range';
    end
    % Orders share a figure only when the saved frequency limit and policy match.
    test = struct('folder',fullfile(root,folders(iFolder).name),'sampleRate',str2double(token.fs),'order',str2double(token.p), ...
        'noisePolicy',noisePolicy,'groupKey',[token.fmax '_' noisePolicy]);
    tests(end+1) = test; %#ok<AGROW>
end
end
