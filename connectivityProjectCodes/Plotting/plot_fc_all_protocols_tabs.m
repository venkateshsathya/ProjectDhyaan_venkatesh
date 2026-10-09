% One figure, one tab per protocol, four occipital seed spectra.
plot_FC = 'granger'; % 'granger' or 'ppc'
[subjects,dates,groups,ages,genders,education,mc] = getDemographicDetails('BK1'); idx = strcmp(subjects,subject);

subject = '054MP';%'012GK';%'040VS';%'096MS';%'019CKa';

disp(table(subjects(idx),groups(idx),ages(idx),genders(idx),education(idx),dates(idx),mc(idx),'VariableNames',{'Subject','Group','Age','Gender','EducationYears','ExperimentDate','MenstrualCycle'}));
protocols = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
if strcmp(plot_FC,'granger')
    % folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-08_17-43-56';
    folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-08_20-16-52';
    folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';

    yLabel = 'Mean outgoing Granger causality';
else
    folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/serverCopies/connectivityProjectCodes/savedData';
    yLabel = 'Mean PPC';
end
seeds = [17 18 48 16]; % Oz, O2, POz, O1 in saved actiCap64_UOL order.
names = {'Oz','O2','POz','O1'};
fig = figure('Color','w','Name',[subject ' | ' plot_FC],'NumberTitle','off');
tabs = uitabgroup(fig);
for p = 1:numel(protocols)
    protocol = protocols{p};
    S = load(fullfile(folder,subject,[protocol '_ep_v8_' plot_FC '.mat']));
    C = (S.connPre + S.connPost)/2; % Same Pre/Post averaging as original script.
    freq = S.freqPost;
    for k = seeds
        C(k,k,:) = NaN; % Exclude self-connections.
    end
    FC = squeeze(mean(C(seeds,:,:),2,'omitnan'));
    tab = uitab(tabs,'Title',protocol);
    ax = axes('Parent',tab);
    plot(ax,freq,FC.','LineWidth',1.5);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,yLabel);
    legend(ax,names,'Location','best'); grid(ax,'on');
    title(ax,[subject ' | ' protocol ' | mean of Pre/Post'],'Interpreter','none');
end
