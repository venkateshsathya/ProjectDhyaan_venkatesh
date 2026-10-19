% GUI copy of plot_fc_all_protocols_tabs.m.
% Saved matrices use C(source, target, frequency). Bad channels contain NaNs.
folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_parametric_2026-10-16_20-59-30';
FC_metric = 'granger_parametric';
fileEnding = ['_ep_v8_' FC_metric '.mat'];
protocols = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
initialSubject = '019CKa';
epochChoice = 'combined'; % 'pre', 'post', or 'combined' = mean of Pre/Post.
frequencyRange = [0 200]; % Hz, displayed in every view.
yLimits = [];            % [] = automatic; or e.g. [0 0.5].

montage = load('actiCap64_UOLLabels.mat');
labels = montage.montageLabels(:,2);
files = dir(fullfile(folder,'*',['*' fileEnding]));
[~,subjects] = cellfun(@fileparts,{files.folder},'UniformOutput',false);
subjects = unique(subjects,'stable');

fig = figure('Color','w','Name','Granger explorer','NumberTitle','off', ...
    'Position',[40 40 1550 900]);
controls = uipanel(fig,'Position',[0 0.92 1 0.08],'BorderType','none');
s.subjectBox = makeDropdown(controls,0.01,'Subject');
s.subjectBox.String = subjects;
if ismember(initialSubject,subjects)
    s.subjectBox.Value = find(strcmp(subjects,initialSubject));
end
s.outBox = makeDropdown(controls,0.21,'Outgoing grid: seed');
s.meanBox = makeDropdown(controls,0.41,'Mean curves: electrode');
s.pairABox = makeDropdown(controls,0.61,'Pair: electrode A');
s.pairBBox = makeDropdown(controls,0.81,'Pair: electrode B');
s.outBox.UserData = 17;  % Prefer Oz when it is a good electrode.
s.meanBox.UserData = 17;
s.pairABox.UserData = 17;
s.pairBBox.UserData = 18; % Prefer O2 for the second electrode.

s.folder = folder; s.fileEnding = fileEnding; s.protocols = protocols;
s.labels = labels; s.epochChoice = epochChoice;
s.frequencyRange = frequencyRange; s.yLimits = yLimits;
s.protocolTabs = uitabgroup(fig,'Position',[0 0 1 0.92]);

fig.UserData = s;
s.subjectBox.Callback = @(~,~) selectSubject(fig);
s.outBox.Callback = @(~,~) drawOutgoing(fig);
s.meanBox.Callback = @(~,~) drawMean(fig);
s.pairABox.Callback = @(~,~) drawPair(fig);
s.pairBBox.Callback = @(~,~) drawPair(fig);
selectSubject(fig);

function box = makeDropdown(parent,left,caption)
uicontrol(parent,'Style','text','String',caption,'Units','normalized', ...
    'Position',[left 0.58 0.18 0.3],'HorizontalAlignment','left');
box = uicontrol(parent,'Style','popupmenu','String',{'Select'}, ...
    'Units','normalized','Position',[left 0.12 0.18 0.42]);
end

function s = makeViews(tab)
% Create the plots the first time this protocol is selected.
s.views = uitabgroup(tab,'Position',[0 0 1 1]);
outTab = uitab(s.views,'Title','64 outgoing plots');
meanTab = uitab(s.views,'Title','Mean outgoing / incoming');
pairTab = uitab(s.views,'Title','Electrode pair');

tiles = tiledlayout(outTab,8,8,'TileSpacing','compact','Padding','compact');
s.outAxes = gobjects(64,1); s.outLines = gobjects(64,1);
for k=1:64
    s.outAxes(k) = nexttile(tiles);
    s.outLines(k) = plot(s.outAxes(k),NaN,NaN,'LineWidth',1);
    s.outAxes(k).FontSize = 7;
    grid(s.outAxes(k),'on');
    if k>56, xlabel(s.outAxes(k),'Hz'); end
    if mod(k-1,8)==0, ylabel(s.outAxes(k),'GC'); end
end
s.meanAxes = axes('Parent',meanTab);
s.meanLines = plot(s.meanAxes,NaN,NaN,'b-',NaN,NaN,'r--','LineWidth',1.5);
xlabel(s.meanAxes,'Frequency (Hz)'); ylabel(s.meanAxes,'Mean Granger causality');
legend(s.meanAxes,{'Outgoing: electrode to others','Incoming: others to electrode'});
grid(s.meanAxes,'on');
s.pairAxes = axes('Parent',pairTab);
s.pairLines = plot(s.pairAxes,NaN,NaN,'b-',NaN,NaN,'r--','LineWidth',1.5);
xlabel(s.pairAxes,'Frequency (Hz)'); ylabel(s.pairAxes,'Granger causality');
grid(s.pairAxes,'on');

end

function selectSubject(fig)
s = fig.UserData;
subject = s.subjectBox.String{s.subjectBox.Value};
s.protocolTabs.SelectionChangedFcn = [];
delete(s.protocolTabs.Children);
files = dir(fullfile(s.folder,subject,['*' s.fileEnding]));
available = erase({files.name},s.fileEnding);
for p=1:numel(s.protocols)
    if ismember(s.protocols{p},available)
        uitab(s.protocolTabs,'Title',s.protocols{p});
    end
end
s.protocolTabs.SelectionChangedFcn = @(~,~) selectProtocol(fig);
selectProtocol(fig);
end

function selectProtocol(fig)
s = fig.UserData;
subject = s.subjectBox.String{s.subjectBox.Value};
protocol = s.protocolTabs.SelectedTab.Title;
D = load(fullfile(s.folder,subject,[protocol s.fileEnding]));
switch s.epochChoice
    case 'pre'
        C = D.connPre; s.freq = D.freqPre;
    case 'post'
        C = D.connPost; s.freq = D.freqPost;
    case 'combined'
        C = (D.connPre + D.connPost)/2; s.freq = D.freqPost;
end
% A bad electrode has no finite incoming or outgoing values in the saved FC.
valid = any(isfinite(C),3);
s.good = find(any(valid,2) & any(valid,1).');
for k=1:64
    C(k,k,:) = NaN; % Exclude self-connections from the plots and means.
end
s.C = C;
tab = s.protocolTabs.SelectedTab;
if isempty(tab.UserData), tab.UserData = makeViews(tab); end
s.plots = tab.UserData;
fig.Name = [subject ' | ' protocol ' | ' s.epochChoice ' | Granger'];
fig.UserData = s;
setElectrodes(s.outBox,s.good,s.labels);
setElectrodes(s.meanBox,s.good,s.labels);
setElectrodes(s.pairABox,s.good,s.labels);
drawOutgoing(fig); drawMean(fig); drawPair(fig);
end

function setElectrodes(box,good,labels)
% Keep the previous electrode if it is still good for this subject/protocol.
previous = box.UserData(box.Value);
box.Value = 1;
items = cell(numel(good),1);
for k=1:numel(good)
    items{k} = sprintf('%02d: %s',good(k),labels{good(k)});
end
box.String = items;
box.UserData = good;
selected = find(good==previous,1);
if ~isempty(selected), box.Value = selected; end
end

function drawOutgoing(fig)
s = fig.UserData;
seed = s.outBox.UserData(s.outBox.Value);
for k=1:64
    s.plots.outLines(k).XData = s.freq;
    s.plots.outLines(k).YData = reshape(s.C(seed,k,:),1,[]);
    caption = sprintf('%s -> %s (%d)',s.labels{seed},s.labels{k},k);
    if ~ismember(k,s.good), caption = [caption ' | excluded']; end
    if k==seed, caption = [caption ' | self']; end
    title(s.plots.outAxes(k),caption,'Interpreter','none','FontSize',8);
end
limits = s.yLimits;
if isempty(limits)
    inRange = s.freq>=s.frequencyRange(1) & s.freq<=s.frequencyRange(2);
    values = s.C(seed,s.good,inRange);
    limits = [0 max(1.05*max([0;values(:)],[],'omitnan'),eps)];
end
set(s.plots.outAxes,'XLim',s.frequencyRange,'YLim',limits); % Same scale across all 64.
end

function drawMean(fig)
s = fig.UserData;
electrode = s.meanBox.UserData(s.meanBox.Value);
outgoing = mean(s.C(electrode,s.good,:),2,'omitnan');
incoming = mean(s.C(s.good,electrode,:),1,'omitnan');
set(s.plots.meanLines(1),'XData',s.freq,'YData',reshape(outgoing,1,[]));
set(s.plots.meanLines(2),'XData',s.freq,'YData',reshape(incoming,1,[]));
title(s.plots.meanAxes,[s.labels{electrode} ' | mean over other good electrodes']);
xlim(s.plots.meanAxes,s.frequencyRange);
if isempty(s.yLimits), ylim(s.plots.meanAxes,'auto'); else, ylim(s.plots.meanAxes,s.yLimits); end
end

function drawPair(fig)
s = fig.UserData;
a = s.pairABox.UserData(s.pairABox.Value);
setElectrodes(s.pairBBox,s.good(s.good~=a),s.labels); % Distinct good electrodes.
b = s.pairBBox.UserData(s.pairBBox.Value);
set(s.plots.pairLines(1),'XData',s.freq,'YData',reshape(s.C(a,b,:),1,[]));
set(s.plots.pairLines(2),'XData',s.freq,'YData',reshape(s.C(b,a,:),1,[]));
legend(s.plots.pairAxes,{[s.labels{a} ' -> ' s.labels{b}], ...
    [s.labels{b} ' -> ' s.labels{a}]},'Location','best');
title(s.plots.pairAxes,[s.labels{a} ' and ' s.labels{b} ' | both directions']);
xlim(s.plots.pairAxes,s.frequencyRange);
if isempty(s.yLimits), ylim(s.plots.pairAxes,'auto'); else, ylim(s.plots.pairAxes,s.yLimits); end
end
