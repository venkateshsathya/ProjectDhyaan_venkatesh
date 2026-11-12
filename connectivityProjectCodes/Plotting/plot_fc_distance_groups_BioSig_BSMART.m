% Overlay all saved BioSig estimators and BSMART, with BioSig 8 excluded.
% Same averaging/group modes as plot_fc_distance_groups_Allprotols_seedavg_compare_runs:
% target mean per distance bin -> seed mean per subject -> equal-weight subject mean.
% Run this script with the MATLAB Run button. No connectivity is recomputed.
%
% BioSig 1/6: multichannel Yule-Walker (biased/unbiased).
% BioSig 2/5: Vieira-Morf (unbiased/biased). BioSig 7: Nuttall-Strand.
% PPC uses its own right y-axis and saved frequency grid; no rescaling of values.
% BSMART is black dash-dot when PPC is enabled; PPC is black dashed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared paper-derived rules apply to EVERY selected protocol and figure.
%{
| Setting | Meaning and reason |
|---|---|
| `usePairedSubjects = true` | Starts with the predefined 30 meditator–control pairs. If either member fails eligibility, both are removed, preserving the demographic matching. |
| `matchSubjectsAcrossRuns = true` | Within each protocol, all available GC estimators and PPC use the **same retained participants**. This prevents a method comparison from accidentally becoming a comparison of different subject groups. |
| `minimumGoodTrials = 30` | Requires at least 30 good trials for eligibility. This checks the available data; equalizing trial counts is a separate computation step. |
| `minimumGoodSeeds = 2` | Requires at least two valid electrodes from the selected seed cluster. Bad electrodes remain excluded through NaN-aware averaging. |
| `seed_electrode_group = 'BackOccipital'` | Uses **Oz, O2, POz and O1**, matching the main paper figures. Changing this setting changes the seeds consistently across every figure in this script. |
| `epochChoice = 'paper'` | For **EO1, EC1, M1, EO2 and EC2**, averages the saved Pre and Post connectivity estimates. For **G1, G2 and M2**, uses Post—the approximately 0.25–1.25 s stimulus window. It averages connectivity estimates rather than concatenating the raw windows. |
| Averaging order | First averages target electrodes within each distance bin, then valid seeds within each subject, then subjects with equal weight. Subjects with more trials do not receive greater plotting weight. |
| Original PPC files | Uses the authors' saved PPC results, preserving their computation. This supports comparison with the paper's PPC curves. |
| `includePPCSelfPairs = true` | Reproduces the [authors' nearest-bin averaging](https://github.com/supratimray/ProjectDhyaanBK1Programs/blob/master/connectivityProjectCodes/displayConnDataAllSubjects.m). **PPC of an electrode with itself equals 1**, so including it raises the closest-bin average. This is a paper-reproduction choice; GC always excludes self-connections. |
| Frequency limits | PPC retains its saved **0–200 Hz** range. Your GC data use **Fs = 200 Hz**, so GC ends at **100 Hz**. PPC's frequency range should not be transferred to these GC results. |
| `ppcOmitLineNoise = true` | Shows gaps within ±3 Hz of 50 Hz and its harmonics. This changes the display; it does not filter or recompute PPC. |
| Automatic GC ymax | Uses **3–95 Hz**, including in the new BioSig 2 comparison. Full saved curves remain plotted. Negative differences remain supported, and explicit `yLimits` override automatic scaling. |
%}

usePairedSubjects = true;               % Keep both members of each valid predefined matched pair.
matchSubjectsAcrossRuns = true;        % Same cohort across the available GC estimators and PPC in each protocol.
minimumGoodTrials = 30;                 % Eligibility cutoff, not a trial-count equalization step.
minimumGoodSeeds = 2;                   % At least two valid seeds from the selected cluster.
includePPCSelfPairs = true;             % Reproduce the authors' nearest-bin PPC average; GC always excludes self-pairs.
ppcOmitLineNoise = true;                % Display-only gaps: +/-3 Hz around 50 Hz and its harmonics.
plotPPCGroupsTogether = true;           % PPC group comparison, with one protocol tab at a time.
plotBioSig2GroupsTogether = true;       % Additional BioSig 2 GC group comparison in protocol tabs.
biosigRoot = '/home/venkatesh/Documents/IISc/data/savedata_parametric_fast_biosig12_2026-11-05_18-37-24';
bsmartFolder = '/home/venkatesh/Documents/IISc/data/100HzMaxFreq_Bsmart_RepeatRunByMistake';
overlayPPC = true;
ppcDataRoot = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/savedData'; %Data from Vinod
% Alternative preprocessing: '/home/venkatesh/Documents/IISc/data/PPC_LineNoiseRemoved'.
ppcYLimits = [];                        % [] = automatic per panel; or [0 1].
biosigOptions = [2];                     % [] = all matching options except 8; or [1 2 5 6 7].
sampleRate = 200;
mvarOrder = 20;
maxFrequency = 100;
ppcMaxFrequency = 100;                  % PPC uses its native saved grid; GC remains limited to maxFrequency.
protocols = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};                     % Or {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'}.
epochChoice = 'paper';                  % Combined for non-stimulus protocols; Post for G1, G2, M2.
% Override with 'pre', 'post', or 'combined' to use the same epoch in every protocol.
subjectsOverride = {};                 % {} = group subjects; or {'053DR'} / a subject list.
plot_modes = {'meditation', 'control', 'combined', 'meditation-control'};              % 'meditation', 'control', 'combined', 'meditation-control'.
% Example: plot_modes = {'meditation','control','combined','meditation-control'};
% Individual subject: set usePairedSubjects=false, subjectsOverride={'053DR'}, plot_modes={'combined'}.
%{
| Dropdown choice | Electrode indices | Electrode names, in the same order |
|---|---|---|
| **BackOccipital** | `[17 18 48 16]` | **Oz, O2, POz, O1** |
| **LeftOccipital** | `[14 44 47]` | **P3, P1, PO3** |
| **RightOccipital** | `[19 49 52]` | **P4, PO4, P2** |
| **Frontal** | `[2 34 35 62 36 63]` | **Fz, AF3, AFz, AF4, F1, F2** |
%}
seed_electrode_dict = dictionary(["BackOccipital", "LeftOccipital","RightOccipital","Frontal"],...
    {[17 18 48 16], [14 44 47],[19 49 52],[2 34 35 62 36 63]});
seed_electrode_group = 'BackOccipital';% ['BackOccipital', 'LeftOccipital', 'RightOccipital', 'Frontal']

badEyeCondition = 'ep';
badTrialVersion = 'v8';
yLimits = [];                          % [] = automatic; or [0 0.5].
autoYMinFrequency = 3;                  % Ignore f < 3 Hz when finding automatic GC ymax.
autoYMaxFrequency = 95;                 % Ignore f > 95 Hz when finding automatic GC ymax.
xLimits = [0 100];                      % Main overlays: PPC to 200 Hz; saved Fs=200-Hz GC to 100 Hz.
gridSetting = 'on';
lineWidth = 1.6;

seeds = seed_electrode_dict{string(seed_electrode_group)};
numSeeds = numel(seeds);
seedAverageRow = numSeeds + 1;

%%%%%%%%%%%%%%%%%%%%%%%% Discover comparison runs %%%%%%%%%%%%%%%%%%%%%%%%
if plotBioSig2GroupsTogether && ~isempty(biosigOptions)
    biosigOptions = unique([biosigOptions 2],'stable');
end
bioTests = findBioSigComparisonRuns(biosigRoot,biosigOptions,sampleRate,mvarOrder,maxFrequency);
folders = {bioTests.folder};
runLabels = arrayfun(@(test) sprintf('BioSig %d',test.algorithm),bioTests,'UniformOutput',false);
folders{end+1} = bsmartFolder;
runLabels{end+1} = 'BSMART';
numRuns = numel(folders);
FC_metric = repmat({'granger_parametric'},1,numRuns);
maxFreqHz = repmat(maxFrequency,1,numRuns);
lineStyles = repmat({'-'},1,numRuns);
lineStyles{end} = '-';
runColors = [lines(numel(bioTests)); 0 0 0];
numGCRuns = numRuns;
if overlayPPC
    lineStyles{end} = '-.';             % Distinguish BSMART from dashed PPC.
    folders{end+1} = ppcDataRoot;
    runLabels{end+1} = 'PPC';
    FC_metric{end+1} = 'ppc';
    maxFreqHz(end+1) = ppcMaxFrequency;
    lineStyles{end+1} = '-';
    runColors(end+1,:) = [0 0 0];
end
fprintf('Comparing Fs=%g Hz, p=%d: %s\n',sampleRate,mvarOrder,strjoin(runLabels,', '));
% Each estimator retains its available subjects and native saved frequency grid.
% Missing/unreadable files are skipped and counts are shown separately.
% Rall.<protocol>.runs{runIndex} retains subject curves, means, counts and skips.
% Rall.<protocol>.ppc retains PPC with the same group/epoch/subject selection.
% figs(m) corresponds to plot_modes{m}; protocol tabs contain six distance bins.

gridSetting = validatestring(gridSetting,{'on','off'});

folders = cellstr(folders); runLabels = cellstr(runLabels); lineStyles = cellstr(lineStyles);
FC_metric = cellstr(FC_metric);
numRuns = numel(folders);
if numRuns < 1 || numel(runLabels) ~= numRuns || numel(lineStyles) ~= numRuns || numel(FC_metric) ~= numRuns
    error('Provide matching folders, runLabels, FC_metric entries and lineStyles.');
end
if any(cellfun(@(metric) isempty(strtrim(metric)),FC_metric))
    error('Each FC_metric entry must be a nonempty saved filename suffix.');
end
if any(~ismember(lineStyles,{'-', '--', ':', '-.'}))
    error('Each lineStyles entry must be ''-'', ''--'', '':'' or ''-.''.');
end
validateattributes(maxFreqHz,{'numeric'},{'vector','numel',numRuns,'real','positive','nonnan'});
validateattributes(lineWidth,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(runColors,{'numeric'},{'size',[numRuns 3],'real','finite','>=',0,'<=',1});
for runIdx = 1:numRuns
    if ~isfolder(folders{runIdx})
        warning('FCOverlay:MissingFolder','%s: folder not found: %s',runLabels{runIdx},folders{runIdx});
    end
end

if ~ismember(lower(epochChoice),{'pre','post','combined','paper'})
    error('epochChoice must be ''pre'', ''post'', ''combined'' or ''paper'' (got ''%s'').',epochChoice);
end
plot_modes = cellstr(plot_modes); plot_modes = plot_modes(:).';
if isempty(plot_modes), error('plot_modes is empty; give at least one mode.'); end
numModes = numel(plot_modes);
modePos = cell(1,numModes); modeNeg = cell(1,numModes);
modeLabel = cell(1,numModes); modeKey = cell(1,numModes);
for m = 1:numModes
    [modePos{m},modeNeg{m},modeLabel{m}] = parseMode(plot_modes{m});
    modeKey{m} = [modePos{m} '-' modeNeg{m}];
end
[~,keepIdx] = unique(modeKey,'stable'); % Ignore repeated modes.
plot_modes = plot_modes(keepIdx); modePos = modePos(keepIdx); modeNeg = modeNeg(keepIdx);
modeLabel = modeLabel(keepIdx); numModes = numel(plot_modes);

[allSubjects,meditationSubjects,controlSubjects] = getGoodSubjectsBK1;
groupSubjects.meditation = meditationSubjects(:);
groupSubjects.control    = controlSubjects(:);
groupSubjects.combined   = allSubjects(:);
pairedSubjects = cell(0,2);
if usePairedSubjects
    pairedSubjects = getPairedSubjectsBK1;
    if ~isempty(subjectsOverride)
        pairedSubjects = pairedSubjects(all(ismember(pairedSubjects,cellstr(subjectsOverride)),2),:);
    end
    groupSubjects.meditation = pairedSubjects(:,1);
    groupSubjects.control = pairedSubjects(:,2);
    groupSubjects.combined = pairedSubjects(:);
end
if ~isempty(subjectsOverride)
    subjectsOverride = cellstr(subjectsOverride);
    groupSubjects.meditation = intersect(groupSubjects.meditation,subjectsOverride,'stable');
    groupSubjects.control = intersect(groupSubjects.control,subjectsOverride,'stable');
    groupSubjects.combined = intersect(groupSubjects.combined,subjectsOverride,'stable');
end

% Groups needed by the requested modes, and the subjects to load for them.
neededGroups = unique([modePos modeNeg(~cellfun(@isempty,modeNeg))],'stable');
if usePairedSubjects || plotPPCGroupsTogether || plotBioSig2GroupsTogether
    neededGroups = unique([neededGroups {'meditation','control'}],'stable');
end
subjectsToLoad = cell(0,1);
for gi = 1:numel(neededGroups)
    subjectsToLoad = [subjectsToLoad; groupSubjects.(neededGroups{gi})]; %#ok<AGROW>
end
subjectsToLoad = unique(subjectsToLoad,'stable');

[electrodeGroups,binNames,binCenters] = getElectrodeGroupsConn('rel',seeds,'actiCap64_UOL');
electrodeGroupsWithSelf = electrodeGroups;
binCenters = fliplr(binCenters); % Match the helper's near-to-far group order.
for r = 1:length(seeds)
    for b = 1:6
        electrodeGroups{r,b} = setdiff(electrodeGroups{r,b},seeds(r),'stable');
    end
end

% One figure per mode, each holding a tab group with one tab per protocol.
figs = gobjects(1,numModes); tabs = cell(numModes,numel(protocols));
for m = 1:numModes
    figs(m) = figure('Color','w','Name',sprintf('BioSig vs BSMART | %s | %s | Seeds: %s | Fs=%g Hz, p=%d',modeLabel{m},epochChoice,seed_electrode_group,sampleRate,mvarOrder), ...
        'NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]);
    tabGroup = uitabgroup(figs(m));
    for p = 1:numel(protocols)
        tabs{m,p} = uitab(tabGroup,'Title',protocols{p});
    end
end
Rall = struct;

for p = 1:numel(protocols)
    protocol = protocols{p};
    thisEpochChoice = epochChoice;
    if strcmpi(epochChoice,'paper')
        if ismember(protocol,{'G1','G2','M2'}), thisEpochChoice = 'post';
        else, thisEpochChoice = 'combined'; end
    end

    runResults = cell(1,numRuns);
    for runIdx = 1:numRuns
        folder = folders{runIdx};
        metric = FC_metric{runIdx};
        runElectrodeGroups = electrodeGroups;
        if strcmp(metric,'ppc') && includePPCSelfPairs
            runElectrodeGroups = electrodeGroupsWithSelf;
        end
        R = struct('folder',folder,'label',runLabels{runIdx},'metric',metric,'lineStyle',lineStyles{runIdx}, ...
            'freq',[],'binNames',{binNames},'binCenters',binCenters, ...
            'skipped',{{}},'subjectNames',{{}},'subjectFC',[], ...
            'epochChoice',thisEpochChoice,'pairedSubjects',{cell(0,2)});
        R.modes = struct('mode',{},'label',{},'meanFC',{},'title',{},'countLabel',{},'message',{});
        for m = 1:numModes
            R.modes(m) = struct('mode',plot_modes{m},'label',modeLabel{m},'meanFC',[], ...
                'title','','countLabel','','message',sprintf('%s: no usable subjects',runLabels{runIdx}));
        end

        % Locate saved subject files for this protocol and run.
        files = dir(fullfile(folder,'*',[protocol '_' badEyeCondition '_' badTrialVersion '_' metric '.mat']));
        availableSubjects = cell(1,numel(files));
        for s = 1:numel(files)
            [~,availableSubjects{s}] = fileparts(files(s).folder);
        end
        [subjectNames,~,inputIndex] = intersect(subjectsToLoad,availableSubjects,'stable');

        % Subjects with no saved file for this protocol.
        skipped = cell(0,2);
        missingFile = setdiff(subjectsToLoad(:),availableSubjects(:),'stable');
        for s = 1:numel(missingFile)
            skipped(end+1,:) = {missingFile{s},'no saved file in folder for this protocol'}; %#ok<SAGROW>
        end

        % Epoch mean -> distance-group target mean -> within-subject seed mean.
        subjectFC = []; freq = [];
        isUsed = false(numel(subjectNames),1);
        for s = 1:numel(subjectNames)
            filePath = fullfile(files(inputIndex(s)).folder,files(inputIndex(s)).name);
            try
                S = load(filePath,'connPre','connPost','freqPre','freqPost','numGoodTrials');
            catch readError
                skipped(end+1,:) = {subjectNames{s},['cannot read file: ' readError.message]}; %#ok<SAGROW>
                continue
            end
            if ~all(isfield(S,{'connPre','connPost','freqPre','freqPost','numGoodTrials'}))
                skipped(end+1,:) = {subjectNames{s},'incomplete saved fields'}; %#ok<SAGROW>
                continue
            end

            % Skip subjects with no usable data (numGoodTrials==0 saves empty arrays).
            if isempty(S.numGoodTrials) || S.numGoodTrials < max(1,minimumGoodTrials) || isempty(S.connPre) || isempty(S.connPost)
                skipped(end+1,:) = {subjectNames{s},sprintf('empty data or fewer than %d good trials',minimumGoodTrials)}; %#ok<SAGROW>
                continue
            end

            switch lower(thisEpochChoice)
                case 'pre', C = S.connPre; thisFreq = S.freqPre;
                case 'post', C = S.connPost; thisFreq = S.freqPost;
                case 'combined'
                    if ~isequal(S.freqPre(:),S.freqPost(:))
                        skipped(end+1,:) = {subjectNames{s},'Pre/Post frequency grids differ'}; %#ok<SAGROW>
                        continue
                    end
                    C = (S.connPre + S.connPost)/2; thisFreq = S.freqPost;
            end

            if isempty(C) || isempty(thisFreq), continue; end
            numGoodSeeds = sum(isfinite(diag(C(seeds,seeds,1))));
            if numGoodSeeds < minimumGoodSeeds
                skipped(end+1,:) = {subjectNames{s},sprintf('%d good seeds; need at least %d',numGoodSeeds,minimumGoodSeeds)}; %#ok<SAGROW>
                continue
            end
            % Allocate on the first VALID subject, using its frequency grid.
            if isempty(subjectFC)
                freq = thisFreq;
                subjectFC = nan(numel(subjectNames),seedAverageRow,6,numel(freq));
            elseif numel(thisFreq) ~= numel(freq) || any(abs(thisFreq(:)-freq(:)) > 1e-9)
                skipped(end+1,:) = {subjectNames{s},'frequency grid differs from the first subject'}; %#ok<SAGROW>
                continue
            end

            for r = 1:length(seeds)
                for b = 1:6
                    subjectFC(s,r,b,:) = reshape(mean(C(seeds(r),runElectrodeGroups{r,b},:),2,'omitnan'),1,1,1,[]);
                end
            end
            subjectFC(s,seedAverageRow,:,:) = mean(subjectFC(s,1:numSeeds,:,:),2,'omitnan');

            % Skip subjects where every seed is bad, so no value exists to average.
            if all(isnan(subjectFC(s,:,:,:)),'all')
                subjectFC(s,:,:,:) = NaN;
                skipped(end+1,:) = {subjectNames{s},'all seed electrodes/targets are NaN (bad electrodes)'}; %#ok<SAGROW>
                continue
            end
            isUsed(s) = true;
        end

        % Keep a matched pair only if both members passed this run's QC.
        if usePairedSubjects
            goodPairs = all(ismember(pairedSubjects,subjectNames(isUsed)),2);
            R.pairedSubjects = pairedSubjects(goodPairs,:);
            keepSubjects = R.pairedSubjects(:);
            partnerMissing = find(isUsed & ~ismember(subjectNames(:),keepSubjects));
            for s = partnerMissing(:).'
                skipped(end+1,:) = {subjectNames{s},'matched partner missing or failed quality criteria'}; %#ok<SAGROW>
            end
            isUsed = isUsed & ismember(subjectNames(:),keepSubjects);
        end

        R.skipped = skipped;
        % Keep the other run available even when this run has no usable data.
        if isempty(subjectFC) || ~any(isUsed)
            warning('FCOverlay:NoData','%s: no usable subjects for protocol %s.',runLabels{runIdx},protocol);
            runResults{runIdx} = R;
            continue
        end

        % Keep only subjects actually used, so names and rows stay aligned.
        subjectFC = subjectFC(isUsed,:,:,:);
        subjectNames = subjectNames(isUsed);

        R.freq = freq(:).';
        R.subjectNames = subjectNames;
        R.subjectFC = subjectFC;
        runResults{runIdx} = R;
    end % runIdx
    if matchSubjectsAcrossRuns
        runResults = alignRunSubjects(runResults);
    end
    % Summarize only after final cohort selection, for all group modes/figures.
    for runIdx = 1:numRuns
        R = summarizeRun(runResults{runIdx},groupSubjects,neededGroups, ...
            plot_modes,modePos,modeNeg,modeLabel,protocol);
        fprintf('%s | %s | %s | %s:',R.label,R.metric,protocol,R.epochChoice);
        for gi = 1:numel(neededGroups)
            gn = neededGroups{gi};
            fprintf(' %s %d/%d;',gn,sum(ismember(R.subjectNames,groupSubjects.(gn))),numel(groupSubjects.(gn)));
        end
        fprintf(' (final used/selected)\n');
        for s = 1:size(R.skipped,1)
            fprintf('  SKIPPED %s: %s\n',R.skipped{s,1},R.skipped{s,2});
        end
        runResults{runIdx} = R;
    end
    Rall.(protocol).runs = runResults(1:numGCRuns);
    if overlayPPC, Rall.(protocol).ppc = runResults{end}; end

    % Overlay all estimators in the same six panels for this protocol/mode.
    for m = 1:numModes
        plotGroupTab(tabs{m,p},runResults,m,binNames,runColors,lineWidth,yLimits,ppcYLimits,gridSetting,maxFreqHz,xLimits,autoYMinFrequency,autoYMaxFrequency,protocol,modeLabel{m},thisEpochChoice,sampleRate,mvarOrder,seed_electrode_group,ppcOmitLineNoise);
    end
end
groupPlotSettings = struct('seedGroup',seed_electrode_group,'maxFrequency',ppcMaxFrequency, ...
    'yLimits',ppcYLimits,'omitLineNoise',ppcOmitLineNoise,'autoYRange',[autoYMinFrequency autoYMaxFrequency], ...
    'gridSetting',gridSetting,'lineWidth',lineWidth,'sampleRate',sampleRate,'mvarOrder',mvarOrder);
if overlayPPC && plotPPCGroupsTogether
    [ppcGroupFigure,ppcGroupTabs] = plotMetricGroupComparison(Rall,protocols,binNames,'ppc',groupPlotSettings);
end
if plotBioSig2GroupsTogether
    groupPlotSettings.maxFrequency = maxFrequency;
    groupPlotSettings.yLimits = yLimits;
    groupPlotSettings.omitLineNoise = false;
    [biosig2GroupFigure,biosig2GroupTabs] = plotMetricGroupComparison(Rall,protocols,binNames,'biosig2',groupPlotSettings);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = alignRunSubjects(results)
% Intersect complete-pair cohorts across runs that have usable data.
% An entirely unavailable run is left empty so the other methods still plot.
available = find(cellfun(@(R) ~isempty(R.subjectNames),results));
available = available(:).';
if isempty(available), return; end
commonSubjects = results{available(1)}.subjectNames;
for k = available(2:end)
    commonSubjects = intersect(commonSubjects,results{k}.subjectNames,'stable');
end
for k = available
    R = results{k};
    keep = ismember(R.subjectNames,commonSubjects);
    for s = find(~keep(:)).'
        R.skipped(end+1,:) = {R.subjectNames{s},'not usable in every available estimator/PPC run'};
    end
    R.subjectNames = R.subjectNames(keep);
    R.subjectFC = R.subjectFC(keep,:,:,:);
    if ~isempty(R.pairedSubjects)
        R.pairedSubjects = R.pairedSubjects(all(ismember(R.pairedSubjects,commonSubjects),2),:);
    end
    results{k} = R;
end
end

function R = summarizeRun(R,groupSubjects,neededGroups,plotModes,modePos,modeNeg,modeLabel,protocol)
if isempty(R.subjectNames), return; end
seedRows = size(R.subjectFC,2);
bins = size(R.subjectFC,3);
nGroup = struct;
for gi = 1:numel(neededGroups)
    gn = neededGroups{gi};
    mask = ismember(R.subjectNames,groupSubjects.(gn));
    nGroup.(gn) = sum(mask);
    R.(gn).subjectNames = R.subjectNames(mask);
    R.(gn).meanFC = reshape(mean(R.subjectFC(mask,:,:,:),1,'omitnan'),seedRows,bins,[]);
    R.(gn).validSubjectCounts = reshape(sum(~isnan(R.subjectFC(mask,:,:,:)),1),seedRows,bins,[]);
end
for m = 1:numel(plotModes)
    pos = modePos{m}; neg = modeNeg{m};
    if isempty(neg)
        R.modes(m).countLabel = sprintf('n=%d',nGroup.(pos));
        if nGroup.(pos)==0, continue; end
        R.modes(m).meanFC = R.(pos).meanFC;
    else
        R.modes(m).countLabel = sprintf('n=%d/%d',nGroup.(pos),nGroup.(neg));
        if nGroup.(pos)==0 || nGroup.(neg)==0, continue; end
        R.modes(m).meanFC = R.(pos).meanFC - R.(neg).meanFC;
    end
    R.modes(m).title = sprintf('%s [%s]: %s | %s | %s | %s', ...
        R.label,R.metric,protocol,modeLabel{m},R.epochChoice,R.modes(m).countLabel);
    R.modes(m).message = '';
end
end

function [pos,neg,label] = parseMode(modeStr)
% 'meditation' -> pos='meditation', neg=''. 'meditation-control' -> pos, neg both set.
s = regexprep(lower(strtrim(char(modeStr))),'\s+','');
parts = strsplit(s,'-');
if numel(parts) > 2 || any(cellfun(@isempty,parts))
    error('Bad plot mode ''%s''. Use ''name'' or ''name-name'' (e.g. ''meditation-control'').',modeStr);
end
pos = canonicalGroup(parts{1},modeStr);
neg = '';
if numel(parts) == 2
    neg = canonicalGroup(parts{2},modeStr);
    if strcmp(pos,neg)
        error('Bad plot mode ''%s'': both sides are the same group.',modeStr);
    end
    label = [displayName(pos) ' - ' displayName(neg)];
else
    label = displayName(pos);
end
end

function g = canonicalGroup(name,modeStr)
switch name
    case {'meditation','meditator','meditators','med'}, g = 'meditation';
    case {'control','controls','con'},                 g = 'control';
    case {'combined','combination','all'},             g = 'combined';
    otherwise
        error('Unknown group ''%s'' in plot mode ''%s''. Use meditation, control or combined.',name,modeStr);
end
end

function d = displayName(g)
switch g
    case 'meditation', d = 'Meditator';
    case 'control',    d = 'Control';
    otherwise,         d = 'Combined';
end
end

function plotGroupTab(tab,runResults,modeIndex,binNames,runColors,lineWidth,yLimits,ppcYLimits,gridSetting,maxFreqHz,xLimits,autoYMinFrequency,autoYMaxFrequency,protocol,modeLabel,epochChoice,sampleRate,mvarOrder,seedGroupLabel,ppcOmitLineNoise)
% Six distance bins; plot only the within-subject seed average for each run.
numRuns = numel(runResults);
hasData = false(1,numRuns);
isPPC = false(1,numRuns);
titleLines = cell(1,numRuns);
for runIdx = 1:numRuns
    mode = runResults{runIdx}.modes(modeIndex);
    hasData(runIdx) = ~isempty(mode.meanFC);
    isPPC(runIdx) = strcmp(runResults{runIdx}.metric,'ppc');
    if hasData(runIdx)
        titleLines{runIdx} = mode.title;
    else
        titleLines{runIdx} = mode.message;
    end
end
if ~any(hasData)
    showMessage(tab,strjoin(titleLines,sprintf('\n')));
    return
end
hasPPC = any(hasData & isPPC);
tl = tiledlayout(tab,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,6);
gcRange = [Inf -Inf];
gcFallbackMax = -Inf;
for b = 1:6
    ax(b) = nexttile(tl); hold(ax(b),'on');
    for runIdx = 1:numRuns
        if ~hasData(runIdx), continue; end
        R = runResults{runIdx};
        meanFC = R.modes(modeIndex).meanFC;
        % The last row was averaged across all selected seeds per subject.
        % Limit each curve on its own saved frequency grid; do not interpolate.
        freqMask = R.freq <= maxFreqHz(runIdx);
        plotFreq = R.freq(freqMask);
        plotValues = reshape(meanFC(end,b,freqMask),1,[]);
        if isPPC(runIdx) && ppcOmitLineNoise, plotValues = maskPPCLineNoise(plotFreq,plotValues); end
        if ~isPPC(runIdx)
            limitMask = isfinite(plotValues);
            if ~isempty(xLimits)
                limitMask = limitMask & plotFreq >= xLimits(1) & plotFreq <= xLimits(2);
            end
            if any(limitMask)
                gcRange(1) = min(gcRange(1),min(plotValues(limitMask)));
                gcFallbackMax = max(gcFallbackMax,max(plotValues(limitMask)));
            end
            upperMask = limitMask & plotFreq >= autoYMinFrequency & plotFreq <= autoYMaxFrequency;
            if any(upperMask), gcRange(2) = max(gcRange(2),max(plotValues(upperMask))); end
        end
        curveLabel = sprintf('%s (%s)',R.label,R.modes(modeIndex).countLabel);
        curveWidth = lineWidth;
        curveTag = 'EstimatorComparisonCurve';
        if isPPC(runIdx)
            yyaxis(ax(b),'right'); hold(ax(b),'on');
            curveLabel = sprintf('PPC (%s, right axis)',R.modes(modeIndex).countLabel);
            curveWidth = max(2,lineWidth);
            curveTag = 'PPCOverlay';
        elseif hasPPC
            yyaxis(ax(b),'left');
        end
        plot(ax(b),plotFreq,plotValues, ...
            'Color',runColors(runIdx,:),'LineStyle',R.lineStyle,'LineWidth',curveWidth,'Marker','none', ...
            'DisplayName',curveLabel, ...
            'Tag',curveTag,'UserData',struct('label',R.label,'metric',R.metric,'folder',R.folder));
    end
    if hasPPC
        yyaxis(ax(b),'right');
        ylabel(ax(b),'Seed-averaged PPC');
        ax(b).YAxis(2).Color = [0 0 0];
        if ~isempty(ppcYLimits), ylim(ax(b),ppcYLimits); end
    end
    if hasPPC, yyaxis(ax(b),'left'); end
    ax(b).YAxis(1).Color = [0 0 0];
    title(ax(b),binNames{b}); grid(ax(b),gridSetting);
    xlabel(ax(b),'Frequency (Hz)'); ylabel(ax(b),'Seed-averaged outgoing GC');
    legend(ax(b),'show','Location','best','Interpreter','none','FontSize',8);
end
% Link frequency only so PPC keeps an independent right-hand scale per bin.
% Preserve the original common GC scale across all six left-hand axes.
linkaxes(ax,'x');
yLimits = chooseGCYLimits(gcRange,gcFallbackMax,yLimits);
for b = 1:6
    if hasPPC, yyaxis(ax(b),'left'); end
    if ~isempty(yLimits), ylim(ax(b),yLimits); end
end
if ~isempty(xLimits), xlim(ax(1),xLimits); end
title(tl,sprintf('%s | %s | %s | Seeds: %s | Fs=%g Hz, p=%d', ...
    protocol,modeLabel,epochChoice,seedGroupLabel,sampleRate,mvarOrder),'Interpreter','none');
end

function values = maskPPCLineNoise(freq,values)
% Display gaps only; the saved spectra and group means remain untouched.
if isempty(freq), return; end
for center = 50:50:max(freq)
    values(:,abs(freq-center)<=3) = NaN;
end
end

function yLimits = chooseGCYLimits(dataRange,fallbackMax,yLimits)
% Only the upper limit ignores the configured frequency edges.
if ~isempty(yLimits), return; end
if ~isfinite(dataRange(2)), dataRange(2) = fallbackMax; end
if ~all(isfinite(dataRange)), return; end
yLimits = [min(0,dataRange(1)) max(0,dataRange(2))];
span = diff(yLimits);
if span==0
    yLimits = [0 1];
else
    if yLimits(1)<0, yLimits(1) = yLimits(1)-0.05*span; end
    yLimits(2) = yLimits(2)+0.05*span;
end
end

function [fig,tabs] = plotMetricGroupComparison(Rall,protocols,binNames,metric,settings)
% Magenta meditators and blue controls share one y-axis in each panel.
% Protocol tabs keep all eight protocols readable in both comparison figures.
isPPC = strcmp(metric,'ppc');
if isPPC
    metricLabel = 'PPC'; curveTag = 'PPCGroupMean'; yLabel = 'Seed-averaged PPC';
else
    metricLabel = sprintf('BioSig 2 GC | Fs=%g Hz, p=%d',settings.sampleRate,settings.mvarOrder);
    curveTag = 'BioSig2GroupMean'; yLabel = 'Seed-averaged outgoing GC';
end
fig = figure('Color','w','Name',[metricLabel ': Meditator vs Control | Seeds: ' settings.seedGroup], ...
    'NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]);
tabGroup = uitabgroup(fig);
tabs = gobjects(1,numel(protocols));
groups = {'meditation','control'};
groupLabels = {'Meditator','Control'};
colors = [0.8 0 0.8;0.25 0.41 0.88];
for p = 1:numel(protocols)
    protocol = protocols{p};
    tabs(p) = uitab(tabGroup,'Title',protocol);
    R = struct;
    if isPPC
        R = Rall.(protocol).ppc;
    else
        candidates = Rall.(protocol).runs;
        index = find(cellfun(@(run) strcmp(run.label,'BioSig 2'),candidates),1);
        if ~isempty(index), R = candidates{index}; end
    end
    if ~isfield(R,'epochChoice') || isempty(R.subjectNames)
        showMessage(tabs(p),sprintf('No usable %s data for %s.',metricLabel,protocol));
        continue
    end
    tl = tiledlayout(tabs(p),3,2,'TileSpacing','compact','Padding','compact');
    axesList = gobjects(1,6);
    gcRange = [Inf -Inf]; fallbackMax = -Inf;
    for b = 1:6
        ax = nexttile(tl); axesList(b) = ax; hold(ax,'on');
        for g = 1:2
            group = groups{g};
            if ~isfield(R,group) || isempty(R.(group).subjectNames), continue; end
            freqMask = R.freq <= settings.maxFrequency;
            freq = R.freq(freqMask);
            values = reshape(R.(group).meanFC(end,b,freqMask),1,[]);
            if isPPC && settings.omitLineNoise, values = maskPPCLineNoise(freq,values); end
            if ~isPPC
                finite = isfinite(values);
                if any(finite)
                    gcRange(1) = min(gcRange(1),min(values(finite)));
                    fallbackMax = max(fallbackMax,max(values(finite)));
                end
                upper = finite & freq>=settings.autoYRange(1) & freq<=settings.autoYRange(2);
                if any(upper), gcRange(2) = max(gcRange(2),max(values(upper))); end
            end
            plot(ax,freq,values,'Color',colors(g,:),'LineWidth',settings.lineWidth,'Marker','none', ...
                'DisplayName',sprintf('%s (n=%d)',groupLabels{g},numel(R.(group).subjectNames)), ...
                'Tag',curveTag,'UserData',struct('group',group,'protocol',protocol));
        end
        title(ax,binNames{b});
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,yLabel); grid(ax,settings.gridSetting);
        xlim(ax,[0 settings.maxFrequency]);
        if ~isempty(ax.Children), legend(ax,'show','Location','best','FontSize',8); end
    end
    limits = settings.yLimits;
    if ~isPPC, limits = chooseGCYLimits(gcRange,fallbackMax,limits); end
    if ~isempty(limits), set(axesList,'YLim',limits); end
    linkaxes(axesList,'x');
    title(tl,sprintf('%s | %s | %s | Seeds: %s',metricLabel,protocol,R.epochChoice,settings.seedGroup),'Interpreter','none');
end
end

function showMessage(tab,msg)
% Empty-tab placeholder when no estimator has data for this protocol/mode.
ax = axes('Parent',tab,'Visible','off','Position',[0 0 1 1]);
text(ax,0.5,0.5,msg,'Units','normalized','HorizontalAlignment','center', ...
    'FontSize',12,'Interpreter','none');
end

function tests = findBioSigComparisonRuns(root,selectedOptions,sampleRate,order,maxFrequency)
tests = struct('folder',{},'algorithm',{});
entries = dir(root);
for k=1:numel(entries)
    if ~entries(k).isdir, continue; end
    token = regexp(entries(k).name, ...
        '^biosig(?<algorithm>\d+)_fs(?<fs>[0-9.eE+-]+)Hz_p(?<p>\d+)_fmax(?<fmax>[0-9.eE+-]+)Hz_withoutLineNoise_keeptrials_no$', ...
        'names','once');
    if isempty(token), continue; end
    algorithm = str2double(token.algorithm);
    if algorithm==8 || (~isempty(selectedOptions) && ~ismember(algorithm,selectedOptions)), continue; end
    if str2double(token.fs)~=sampleRate || str2double(token.p)~=order || str2double(token.fmax)~=maxFrequency
        continue
    end
    tests(end+1) = struct('folder',fullfile(root,entries(k).name),'algorithm',algorithm); %#ok<AGROW>
end
[~,index] = sort([tests.algorithm]);
tests = tests(index);
end
