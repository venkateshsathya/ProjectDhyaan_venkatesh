% Six distance bins in a 3-by-2 grid, with all four seeds and Combined overlaid.
% GOOD-SUBJECTS VERSION of plot_fc_distance_groups.m: subjects whose saved file is
% missing or empty (numGoodTrials==0, empty connPre/connPost) are skipped and
% reported in the Command Window, instead of crashing the script.
% Edit these inputs, then run this script with the MATLAB Run button.
folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';
protocol = 'EO2';          % EO1, EC1, G1, M1, G2, EO2, EC2, or M2.
subjectGroup = 'meditation'; % 'meditation', 'control', or 'combined'.
epochChoice = 'post';  % 'pre', 'post', or 'combined'.

% DATA AND SUBJECT SELECTION
% Saved matrices are C(source,target,frequency), in actiCap64_UOL order.
% All supplied subjects and requested epochs must share the frequency grid.
% Thus C(seed,j,f) is outgoing directed GC from the seed to target j.
% No symmetrization, absolute value, or restriction to [0,1] is applied.
% getGoodSubjectsBK1 removes the project's declared bad subjects.
% Only selected subjects with saved, non-empty files in folder are used.
% Subjects that are selected but skipped are listed in R.skipped
% (name + reason) and printed to the Command Window.
% Saved NaNs remove bad electrodes. Self-connections are also excluded.
% No additional trial-count cutoff or negative-value rejection is applied.
% The project's information files and montage must be on the MATLAB path.
%
% DISTANCE GROUPS (SAME DEFINITION AS getElectrodeGroupsConn)
% delta_a(r,j), delta_e(r,j): shortest azimuth/elevation differences in degrees.
% theta(r,j) = sqrt(delta_a(r,j)^2 + delta_e(r,j)^2).
% x(r,j) = cos(pi*theta(r,j)/180); the helper sets delta_a(r,Cz)=0.
% Six equal-width bins cover [-1,1]; larger x means nearer electrodes.
% Each seed has its OWN target membership B(r,b), excluding itself.
% The helper returns groups from near to far, so centers are reversed below.
% This is the project's angular-space definition, not a 3-D dot product.
%
% AVERAGING MATH (ALL SUMS ARE AT EACH FREQUENCY f)
% s = subject, r = seed, j = target electrode, b = distance bin.
% E = requested epochs: {pre}, {post}, or {pre,post}.
% X(s,r,j,f) = (1/|E|) sum_{t in E} C(s,r,j,f,t).
% Combined epochs require both pair estimates; a missing one leaves NaN.
% V(s,r,b,f) = {j in B(r,b) : X(s,r,j,f) is not NaN}.
% A(s,r,b,f) = (1/|V|) sum_{j in V(s,r,b,f)} X(s,r,j,f).
% K(s,b,f) = {r in {Oz,O2,POz,O1} : A(s,r,b,f) is not NaN}.
% A(s,combined,b,f) = (1/|K|) sum_{r in K(s,b,f)} A(s,r,b,f).
% U(r,b,f) = {s in the selected subject group : A(s,r,b,f) is not NaN}.
% Y(r,b,f) = (1/|U|) sum_{s in U(r,b,f)} A(s,r,b,f).
% Empty averages remain NaN. Averages weight targets, seeds, and subjects
% equally at their respective stages; subjects with more good electrodes
% receive no extra weight. Combined subjects pool individuals from both
% groups, rather than taking a 50:50 average of the two group means.
% Combined seeds use the available seeds for each subject/bin/frequency.
%
% RESULTS LEFT IN THE WORKSPACE (F=201 for the inspected saved files)
% R.subjectFC         : Nsubjects-by-5-by-6-by-F (only subjects actually used).
% R.meanFC            : 5-by-6-by-F (the plotted population means Y).
% R.validSubjectCounts: 5-by-6-by-F (the denominators |U|).
% R.subjectNames      : names matching the rows of R.subjectFC.
% R.skipped           : Nskipped-by-2 cell, {subjectName, reason}.
% Seed order: Oz, O2, POz, O1, Combined. Bin order: near to far.
% R.electrodeGroups   : 4-by-6 cell array of target-electrode indices.
% R.freq              : 1-by-F; R.binCenters: 1-by-6.
% All five curves are calculated and plotted in every panel.

seeds = [17 18 48 16]; % Oz, O2, POz, O1 in the saved 64-channel order.
seedNames = {'Oz','O2','POz','O1','Combined'};
[allSubjects,meditationSubjects,controlSubjects] = getGoodSubjectsBK1;
switch lower(subjectGroup)
    case 'meditation', selectedSubjects = meditationSubjects;
    case 'control', selectedSubjects = controlSubjects;
    case 'combined', selectedSubjects = allSubjects;
end
[electrodeGroups,binNames,binCenters] = getElectrodeGroupsConn('rel',seeds,'actiCap64_UOL');
binCenters = fliplr(binCenters); % Match the helper's near-to-far group order.
for r = 1:4
    for b = 1:6
        electrodeGroups{r,b} = setdiff(electrodeGroups{r,b},seeds(r),'stable');
    end
end

% Locate saved subject files for the selected protocol.
files = dir(fullfile(folder,'*',[protocol '_ep_v8_granger.mat']));
availableSubjects = cell(1,numel(files));
for s = 1:numel(files)
    [~,availableSubjects{s}] = fileparts(files(s).folder);
end
[subjectNames,~,inputIndex] = intersect(selectedSubjects,availableSubjects,'stable');

% Selected subjects with no saved file for this protocol (e.g. never processed,
% or failed in the saving step). They are recorded as skipped.
skipped = cell(0,2);
missingFile = setdiff(selectedSubjects(:),availableSubjects(:),'stable');
for s = 1:numel(missingFile)
    skipped(end+1,:) = {missingFile{s},'no saved file in folder for this protocol'}; %#ok<SAGROW>
end

% Epoch mean -> distance-group target mean -> within-subject seed mean.
subjectFC = []; freq = []; % Start fresh on each script run.
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

if isempty(subjectFC) || ~any(isUsed)
    error('No subjects with usable data for protocol %s, group %s.',protocol,subjectGroup);
end

% Keep only subjects actually used, so names and rows stay aligned.
subjectFC = subjectFC(isUsed,:,:,:);
subjectNames = subjectNames(isUsed);

% Report what was dropped.
fprintf('%s | %s | %s: %d of %d selected subjects used.\n', ...
    protocol,subjectGroup,epochChoice,numel(subjectNames),numel(selectedSubjects));
for s = 1:size(skipped,1)
    fprintf('  SKIPPED %s: %s\n',skipped{s,1},skipped{s,2});
end

% Subject mean is the LAST averaging step, preserving equal subject weights.
R = struct;
R.protocol = protocol; R.subjectGroup = subjectGroup; R.epochChoice = epochChoice;
R.subjectNames = subjectNames; R.subjectFC = subjectFC; R.freq = freq(:).';
R.skipped = skipped;
R.meanFC = reshape(mean(subjectFC,1,'omitnan'),5,6,[]);
R.validSubjectCounts = reshape(sum(~isnan(subjectFC),1),5,6,[]);
R.seedNames = seedNames; R.seedIndices = seeds;
R.electrodeGroups = electrodeGroups; R.binNames = binNames; R.binCenters = binCenters;

% One panel per distance bin; automatic limits retain the actual GC range.
fig = figure('Color','w','Name',[char(protocol) ' | FC by distance'],'NumberTitle','off');
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,6); colors = [lines(4); 0 0 0];
for b = 1:6
    ax(b) = nexttile; hold(ax(b),'on');
    for r = 1:5
        plot(ax(b),R.freq,reshape(R.meanFC(r,b,:),1,[]),'Color',colors(r,:),'LineWidth',1.5);
    end
    title(ax(b),binNames{b}); grid(ax(b),'on');
    xlabel(ax(b),'Frequency (Hz)'); ylabel(ax(b),'Mean outgoing Granger causality');
    legend(ax(b),seedNames,'Location','best');
end
linkaxes(ax,'xy');
sgtitle(sprintf('%s | %s | %s | N=%d used of %d (meditation=%d, control=%d)', ...
    char(protocol),char(subjectGroup),char(epochChoice),numel(subjectNames),numel(selectedSubjects), ...
    sum(ismember(subjectNames,meditationSubjects)),sum(ismember(subjectNames,controlSubjects))));