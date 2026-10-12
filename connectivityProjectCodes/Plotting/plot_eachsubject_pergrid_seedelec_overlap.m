% Run this script: one subject at a time in the same 4-by-2 figure.
% Each protocol panel overlays Oz, O2, POz, O1. FC is averaged over targets.
folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';
protocols = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
epochChoice = 'combined'; % 'pre', 'post', or 'combined'.
saveVideo = true;       % Set true to record while the figure updates.
% videoFile = fullfile(fileparts(mfilename('fullpath')),'fc_subject_movie_protocols.avi');
timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
videoFile = fullfile(fileparts(mfilename('fullpath')),['fc_subject_movie' timestamp '.avi']);

% Use good subjects with a folder here, or replace subjects with your own list.
subjects = getGoodSubjectsBK1;
D = dir(folder);
subjects = intersect(subjects,{D([D.isdir]).name},'stable');
% subjects = {'019CKa','040VS','035SS'};
seeds = [17 18 48 16]; % Oz, O2, POz, O1 in the saved channel order.
names = {'Oz','O2','POz','O1'};
colors = lines(4);
fig = figure('Color','w','Position',[100 50 1200 900],'Resize','off');

% Get number of hours meditated details
[sn,~,group,age,gender] = getDemographicDetails('BK1');
P = load('BK1PracticeHours.mat');
hours = nan(size(sn));
[found,k] = ismember(sn,P.practiceHours(2:end,1));
hours(found) = cell2mat(P.practiceHours(k(found)+1,2));

% Motion JPEG AVI works with MATLAB on Linux, plays with VLC on a Mac,
% and can be uploaded to YouTube. One recorded frame lasts one second.
if saveVideo
    v = VideoWriter(videoFile,'Motion JPEG AVI');
    v.FrameRate = 1;
    v.Quality = 95;
    open(v);
end

for s = 1:numel(subjects)
    % Load all protocols before drawing this subject's frame.
    FC = cell(1,numel(protocols)); freq = cell(1,numel(protocols));
    for p = 1:numel(protocols)
        S = load(fullfile(folder,subjects{s},[protocols{p} '_ep_v8_granger.mat']));
        FC{p} = nan(4,0); % A zero-trial protocol has no curve.
        if S.numGoodTrials > 0
            switch epochChoice
                case 'pre', C = S.connPre; freq{p} = S.freqPre;
                case 'post', C = S.connPost; freq{p} = S.freqPost;
                case 'combined', C = (S.connPre + S.connPost)/2; freq{p} = S.freqPost;
            end
            for k = seeds
                C(k,k,:) = NaN; % Exclude each seed's self-connection.
            end
            % Outgoing GC: mean over targets; saved bad-channel NaNs are omitted.
            FC{p} = reshape(mean(C(seeds,:,:),2,'omitnan'),4,[]);
        end
    end

    figure(fig); clf(fig);
    for p = 1:numel(protocols)
        subplot(4,2,p); hold on;
        for r = 1:4
            plot(freq{p},FC{p}(r,:),'Color',colors(r,:),'LineWidth',1.2);
        end
        title(protocols{p}); grid on;
        xlabel('Frequency (Hz)'); ylabel('Mean outgoing Granger causality');
        legend(names,'Location','best');
    end
    % sgtitle(sprintf('%s | %s | subject %d of %d',subjects{s},epochChoice,s,numel(subjects)));
    % sgtitle(sprintf('%s | %s | Group: %s | Age: %g | Gender: %s',subjects{s},epochChoice,group{idx},age(idx),gender{idx}));
    % sgtitle(sprintf('%s | %s | subject %d of %d',subjects{s},epochChoice,s,numel(subjects)));
    % [sn,~,group,age,gender] = getDemographicDetails('BK1'); idx = strcmp(sn,subjects{s}); 
    idx = strcmp(sn,subjects{s});
    sgtitle(sprintf('%s | %s | Group: %s | Age: %g | Gender: %s | Meditation: %g hours',subjects{s},epochChoice,group{idx},age(idx),gender{idx},hours(idx)));

    set(findall(fig,'Type','axes'),'YLim',[0 1.5]);
    drawnow;
    if saveVideo
        writeVideo(v,getframe(fig)); % Record the displayed frame in this same loop.
    end
    pause(1); % Display this subject for at least one second.
end

if saveVideo
    close(v);
end
