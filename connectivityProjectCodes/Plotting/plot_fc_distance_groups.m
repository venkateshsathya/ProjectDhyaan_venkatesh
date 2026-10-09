function [R,fig] = plot_fc_distance_groups(protocol,subjectGroup,epochChoice,fcInput,seedChoice)
% PLOT_FC_DISTANCE_GROUPS: six distance bins in a 3-by-2 frequency-plot grid.
%
% INPUTS
% protocol     : 'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', or 'M2'.
% subjectGroup : 'meditation', 'control', or 'combined'.
% epochChoice  : 'pre', 'post', or 'combined'.
% fcInput      : saved-data folder containing subject/protocol MAT files,
%                OR a cell array of loaded subject structures for ONE protocol.
% seedChoice   : 'all' (default: four seeds + combined), 'combined', or
%                'Oz', 'O2', 'POz', 'O1' (plot only that curve).
%
% EXAMPLE USING SAVED FILES
% folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';
% [R,fig] = plot_fc_distance_groups('EO1','meditation','combined',folder,'all');
%
% EXAMPLE USING MATRICES ALREADY IN MATLAB
% D{1} = struct('subjectName','019CKa','connPre',connPre,'connPost',connPost, ...
%               'freqPre',freqPre,'freqPost',freqPost);
% D{2} = anotherSubjectStructure; % Same fields; each matrix is 64-by-64-by-F.
% [R,fig] = plot_fc_distance_groups('M1','combined','post',D,'all');
% For a single-epoch input only its conn/freq fields are needed.
% For several protocols, call this function once for each protocol.
%
% DATA AND SUBJECT SELECTION
% Saved matrices are C(source,target,frequency), in actiCap64_UOL order.
% All supplied subjects and requested epochs must share the frequency grid.
% Thus C(seed,j,f) is outgoing directed GC from the seed to target j.
% No symmetrization, absolute value, or restriction to [0,1] is applied.
% getGoodSubjectsBK1 removes the project's declared bad subjects.
% Only subjects supplied in fcInput are used; their names are returned in R.
% Saved NaNs remove bad electrodes. Self-connections are also excluded.
% This function adds no trial-count cutoff or negative-value rejection.
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
% OUTPUT DIMENSIONS (F=201 for the inspected saved files)
% R.subjectFC         : Nsubjects-by-5-by-6-by-F.
% R.meanFC            : 5-by-6-by-F (the plotted population means Y).
% R.validSubjectCounts: 5-by-6-by-F (the denominators |U|).
% Seed order: Oz, O2, POz, O1, Combined. Bin order: near to far.
% R.electrodeGroups   : 4-by-6 cell array of target-electrode indices.
% R.freq              : 1-by-F; R.binCenters: 1-by-6.
% Changing seedChoice only selects displayed curves; R keeps all five.

if nargin < 5, seedChoice = 'all'; end
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

% Locate supplied subjects before averaging; saved files are read one at a time.
fromFiles = ischar(fcInput) || isstring(fcInput);
if fromFiles
    files = dir(fullfile(fcInput,'*',[char(protocol) '_ep_v8_granger.mat']));
    availableSubjects = cell(1,numel(files));
    for s = 1:numel(files)
        [~,availableSubjects{s}] = fileparts(files(s).folder);
    end
else
    availableSubjects = cellfun(@(D) D.subjectName,fcInput,'UniformOutput',false);
end
[subjectNames,~,inputIndex] = intersect(selectedSubjects,availableSubjects,'stable');

% Epoch mean -> distance-group target mean -> within-subject seed mean.
for s = 1:numel(subjectNames)
    if fromFiles
        S = load(fullfile(files(inputIndex(s)).folder,files(inputIndex(s)).name));
    else
        S = fcInput{inputIndex(s)};
    end
    switch lower(epochChoice)
        case 'pre', C = S.connPre; freq = S.freqPre;
        case 'post', C = S.connPost; freq = S.freqPost;
        case 'combined', C = (S.connPre + S.connPost)/2; freq = S.freqPost;
    end
    if s == 1, subjectFC = nan(numel(subjectNames),5,6,numel(freq)); end
    for r = 1:4
        for b = 1:6
            subjectFC(s,r,b,:) = reshape(mean(C(seeds(r),electrodeGroups{r,b},:),2,'omitnan'),1,1,1,[]);
        end
    end
    subjectFC(s,5,:,:) = mean(subjectFC(s,1:4,:,:),2,'omitnan');
end

% Subject mean is the LAST averaging step, preserving equal subject weights.
R.protocol = protocol; R.subjectGroup = subjectGroup; R.epochChoice = epochChoice;
R.subjectNames = subjectNames; R.subjectFC = subjectFC; R.freq = freq(:).';
R.meanFC = reshape(mean(subjectFC,1,'omitnan'),5,6,[]);
R.validSubjectCounts = reshape(sum(~isnan(subjectFC),1),5,6,[]);
R.seedNames = seedNames; R.seedIndices = seeds;
R.electrodeGroups = electrodeGroups; R.binNames = binNames; R.binCenters = binCenters;

% One panel per distance bin; automatic limits retain the actual GC range.
if strcmpi(seedChoice,'all'), curves = 1:5;
else, curves = find(strcmpi(seedNames,seedChoice)); end
fig = figure('Color','w','Name',[char(protocol) ' | FC by distance'],'NumberTitle','off');
tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,6); colors = [lines(4); 0 0 0];
for b = 1:6
    ax(b) = nexttile; hold(ax(b),'on');
    for r = curves
        plot(ax(b),R.freq,reshape(R.meanFC(r,b,:),1,[]),'Color',colors(r,:),'LineWidth',1.5);
    end
    title(ax(b),binNames{b}); grid(ax(b),'on');
    xlabel(ax(b),'Frequency (Hz)'); ylabel(ax(b),'Mean outgoing Granger causality');
    legend(ax(b),seedNames(curves),'Location','best');
end
linkaxes(ax,'xy');
sgtitle(sprintf('%s | %s | %s | N=%d (meditation=%d, control=%d)', ...
    char(protocol),char(subjectGroup),char(epochChoice),numel(subjectNames), ...
    sum(ismember(subjectNames,meditationSubjects)),sum(ismember(subjectNames,controlSubjects))));
end
