% One FIGURE per entry in plot_modes. Each figure has one tab per protocol; each tab
% shows the six distance bins (3-by-2 grid) with all four seeds and Combined overlaid.
% Edit these inputs, then run this script with the MATLAB Run button.

folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';
folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-08_20-16-52';

protocols = {'EO1'};%{'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
epochChoice = 'combined';  % 'pre', 'post', or 'combined'.

% PLOT MODES (a cell array; one new figure per entry, in the order given).
%   'meditation'  -> mean over meditation subjects
%   'control'     -> mean over control subjects
%   'combined'    -> mean over all subjects pooled (individuals, not a 50:50 mean of groups)
%   'A-B'         -> difference of group means, A minus B, for A,B in the three names above
%                    (e.g. 'meditation-control' or 'control-meditation')
% Also accepted: 'controls', 'meditators', 'all'; case and spaces are ignored.
% Examples: plot_modes = {'meditation','control'};
%           plot_modes = {'meditation','control','combined','meditation-control'};
plot_modes = {'control'};%{'meditation','control','combined','meditation-control'};

% CALCULATION (same as plot_fc_distance_groups_claude.m, per protocol)
% Per subject: epoch mean -> target mean per distance bin -> seed mean, giving
% subjectFC(s,seed,bin,f) with seeds Oz, O2, POz, O1, Combined. This does not
% depend on the group, so it is computed ONCE per protocol, only for the subjects
% needed by the requested modes. Group means are then taken over subjects with NaNs
% omitted. A difference is NaN where either group mean is missing and carries no
% statistics (unpaired group means).
% Subjects with a missing or empty saved file (numGoodTrials==0, empty
% connPre/connPost) are skipped for that protocol and listed in the Command Window.
% No abs, no negative-value rejection and no trial/seed cutoff is applied.
%
% RESULTS LEFT IN THE WORKSPACE
% Rall.<protocol>.<group> for each group used by the modes (meditation, control,
%   combined): .meanFC (5-by-6-by-F), .validSubjectCounts (5-by-6-by-F), .subjectNames.
% Rall.<protocol>.modes(m): .mode, .label, .meanFC (5-by-6-by-F, or [] if unavailable)
%   for the m-th entry of plot_modes (differences included).
% Rall.<protocol>.freq, .binNames, .binCenters, .skipped, .subjectNames, .subjectFC.
% figs(m) is the figure handle for plot_modes{m}.

seeds = [17 18 48 16]; % Oz, O2, POz, O1 in the saved 64-channel order.
seedNames = {'Oz','O2','POz','O1','Combined'};

% Validate inputs before doing any work.
if ~ismember(lower(epochChoice),{'pre','post','combined'})
    error('epochChoice must be ''pre'', ''post'' or ''combined'' (got ''%s'').',epochChoice);
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

% Groups needed by the requested modes, and the subjects to load for them.
neededGroups = unique([modePos modeNeg(~cellfun(@isempty,modeNeg))],'stable');
subjectsToLoad = cell(0,1);
for gi = 1:numel(neededGroups)
    subjectsToLoad = [subjectsToLoad; groupSubjects.(neededGroups{gi})]; %#ok<AGROW>
end
subjectsToLoad = unique(subjectsToLoad,'stable');

[electrodeGroups,binNames,binCenters] = getElectrodeGroupsConn('rel',seeds,'actiCap64_UOL');
binCenters = fliplr(binCenters); % Match the helper's near-to-far group order.
for r = 1:4
    for b = 1:6
        electrodeGroups{r,b} = setdiff(electrodeGroups{r,b},seeds(r),'stable');
    end
end

% One figure per mode, each holding a tab group with one tab per protocol.
figs = gobjects(1,numModes); tabs = cell(numModes,numel(protocols));
for m = 1:numModes
    figs(m) = figure('Color','w','Name',sprintf('FC by distance | %s | %s',modeLabel{m},epochChoice), ...
        'NumberTitle','off','Units','normalized','OuterPosition',[0 0 1 1]);
    tabGroup = uitabgroup(figs(m));
    for p = 1:numel(protocols)
        tabs{m,p} = uitab(tabGroup,'Title',protocols{p});
    end
end
colors = [lines(4); 0 0 0];
Rall = struct;

for p = 1:numel(protocols)
    protocol = protocols{p};

    % Locate saved subject files for this protocol.
    files = dir(fullfile(folder,'*',[protocol '_ep_v8_granger.mat']));
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
        S = load(fullfile(files(inputIndex(s)).folder,files(inputIndex(s)).name));

        % Skip subjects with no usable data (numGoodTrials==0 saves empty arrays).
        if S.numGoodTrials == 0 || isempty(S.connPre) || isempty(S.connPost)
            skipped(end+1,:) = {subjectNames{s},'saved file has no data (numGoodTrials==0 or empty conn)'}; %#ok<SAGROW>
            continue
        end

        switch lower(epochChoice)
            case 'pre', C = S.connPre; thisFreq = S.freqPre;
            case 'post', C = S.connPost; thisFreq = S.freqPost;
            case 'combined', C = (S.connPre + S.connPost)/2; thisFreq = S.freqPost;
        end

        % Allocate on the first VALID subject, using its frequency grid.
        if isempty(subjectFC)
            freq = thisFreq;
            subjectFC = nan(numel(subjectNames),5,6,numel(freq));
        elseif numel(thisFreq) ~= numel(freq) || any(abs(thisFreq(:)-freq(:)) > 1e-9)
            skipped(end+1,:) = {subjectNames{s},'frequency grid differs from the first subject'}; %#ok<SAGROW>
            continue
        end

        for r = 1:4
            for b = 1:6
                subjectFC(s,r,b,:) = reshape(mean(C(seeds(r),electrodeGroups{r,b},:),2,'omitnan'),1,1,1,[]);
            end
        end
        subjectFC(s,5,:,:) = mean(subjectFC(s,1:4,:,:),2,'omitnan');

        % Skip subjects where every seed is bad, so no value exists to average.
        if all(isnan(subjectFC(s,:,:,:)),'all')
            subjectFC(s,:,:,:) = NaN;
            skipped(end+1,:) = {subjectNames{s},'all seed electrodes/targets are NaN (bad electrodes)'}; %#ok<SAGROW>
            continue
        end
        isUsed(s) = true;
    end

    % Report what was used and dropped for this protocol.
    fprintf('%s | %s:',protocol,epochChoice);
    for gi = 1:numel(neededGroups)
        gn = neededGroups{gi};
        fprintf(' %s %d/%d;',gn,sum(isUsed(:) & ismember(subjectNames(:),groupSubjects.(gn))),numel(groupSubjects.(gn)));
    end
    fprintf(' (used/selected)\n');
    for s = 1:size(skipped,1)
        fprintf('  SKIPPED %s: %s\n',skipped{s,1},skipped{s,2});
    end

    % Nothing usable at all: say so in every figure and continue.
    if isempty(subjectFC) || ~any(isUsed)
        warning('No subjects with usable data for protocol %s.',protocol);
        for m = 1:numModes
            showMessage(tabs{m,p},sprintf('%s: no usable subjects',protocol));
        end
        continue
    end

    % Keep only subjects actually used, so names and rows stay aligned.
    subjectFC = subjectFC(isUsed,:,:,:);
    subjectNames = subjectNames(isUsed);

    % Group means for the groups needed by the requested modes.
    R = struct('freq',freq(:).','binNames',{binNames},'binCenters',binCenters, ...
        'skipped',{skipped},'subjectNames',{subjectNames},'subjectFC',subjectFC);
    nGroup = struct;
    for gi = 1:numel(neededGroups)
        gn = neededGroups{gi};
        gMask = ismember(subjectNames(:),groupSubjects.(gn));
        nGroup.(gn) = sum(gMask);
        R.(gn).subjectNames = subjectNames(gMask);
        R.(gn).meanFC = reshape(mean(subjectFC(gMask,:,:,:),1,'omitnan'),5,6,[]);
        R.(gn).validSubjectCounts = reshape(sum(~isnan(subjectFC(gMask,:,:,:)),1),5,6,[]);
    end

    % Mean data, title and any problem message for each requested mode.
    modeMean = cell(1,numModes); modeTitle = cell(1,numModes); modeMsg = cell(1,numModes);
    R.modes = struct('mode',{},'label',{},'meanFC',{});
    for m = 1:numModes
        pos = modePos{m}; neg = modeNeg{m};
        if isempty(neg)
            if nGroup.(pos) == 0
                modeMsg{m} = sprintf('%s: no %s subjects',protocol,lower(displayName(pos)));
            else
                modeMean{m} = R.(pos).meanFC;
                modeTitle{m} = sprintf('%s | %s | %s | N=%d',protocol,modeLabel{m},epochChoice,nGroup.(pos));
            end
        else
            if nGroup.(pos) == 0 || nGroup.(neg) == 0
                modeMsg{m} = sprintf('%s: %s needs subjects in both groups',protocol,modeLabel{m});
            else
                modeMean{m} = R.(pos).meanFC - R.(neg).meanFC;
                modeTitle{m} = sprintf('%s | %s | %s | N(%s)=%d, N(%s)=%d',protocol,modeLabel{m},epochChoice, ...
                    displayName(pos),nGroup.(pos),displayName(neg),nGroup.(neg));
            end
        end
        R.modes(m).mode = plot_modes{m}; R.modes(m).label = modeLabel{m}; R.modes(m).meanFC = modeMean{m};
    end
    Rall.(protocol) = R;

    % Plot this protocol's tab in each figure.
    for m = 1:numModes
        if isempty(modeMean{m})
            showMessage(tabs{m,p},modeMsg{m});
        else
            plotGroupTab(tabs{m,p},modeMean{m},R.freq,binNames,seedNames,colors,modeTitle{m});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    case 'meditation', d = 'Meditation';
    case 'control',    d = 'Control';
    otherwise,         d = 'Combined';
end
end

function plotGroupTab(tab,meanFC,freq,binNames,seedNames,colors,titleStr)
% Six distance bins in a 3-by-2 grid; all seeds plotted in every panel.
tl = tiledlayout(tab,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,6);
for b = 1:6
    ax(b) = nexttile(tl); hold(ax(b),'on');
    for r = 1:5
        plot(ax(b),freq,reshape(meanFC(r,b,:),1,[]),'Color',colors(r,:),'LineWidth',1.5);
    end
    title(ax(b),binNames{b}); grid(ax(b),'on');
    xlabel(ax(b),'Frequency (Hz)'); ylabel(ax(b),'Mean outgoing Granger causality');
    legend(ax(b),seedNames,'Location','best');
end
linkaxes(ax,'xy');
title(tl,titleStr);
end

function showMessage(tab,msg)
% Empty-tab placeholder for protocols/modes without usable subjects.
ax = axes('Parent',tab,'Visible','off','Position',[0 0 1 1]);
text(ax,0.5,0.5,msg,'Units','normalized','HorizontalAlignment','center','FontSize',12);
end