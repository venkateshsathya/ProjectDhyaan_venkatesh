% Four occipital seeds: outgoing GC averaged over all other electrodes.

plot_FC = 'granger';
if strcmp(plot_FC,'granger')
    folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';
                    
else
    folder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/serverCopies/connectivityProjectCodes/savedData';
end

subject = '019CKa';
protocol = 'M1';
if strcmp(plot_FC,'granger')
    S = load(fullfile(folder,subject,[protocol '_ep_v8_granger.mat']));
else
    S = load(fullfile(folder,subject,[protocol '_ep_v8_ppc.mat']));

end
% EO1: average Pre/Post. For stimulus-period GC, use C = S.connPost instead.
C = (S.connPre + S.connPost)/2;
freq = S.freqPost;
seeds = [17 18 48 16]; % Oz, O2, POz, O1 in the saved actiCap64_UOL channel order.
% seeds = [5 6 7 8];
names = {'Oz','O2','POz','O1'};
for k = seeds
    C(k,k,:) = NaN; % Exclude self-connections.
end
FC = squeeze(mean(C(seeds,:,:),2,'omitnan')); % Rows are sources; ignore missing electrodes.

figure('Color','w');
plot(freq,FC.','LineWidth',1.5);
xlabel('Frequency (Hz)'); ylabel('Mean outgoing Granger causality');
legend(names,'Location','best'); grid on;
title([subject ' | ' protocol],'Interpreter','none');
