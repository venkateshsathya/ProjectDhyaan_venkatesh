% Run after your normal project/FieldTrip path setup. Requires Parallel Computing Toolbox.
% Subjects run in parallel; protocols and Pre/Post stay sequential within each subject.
parentFolder='/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1';
ftDataFolder=fullfile(parentFolder,'data','ftData');
subjects=getGoodSubjectsBK1;
% subjects={'019CKa'}; % Use this line for a first single-subject run.
protocolNameList={'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
badEyeCondition='ep'; badTrialVersion='v8'; stRange=[0.25 1.25];
numWorkers=2; % Increase after checking runtime and RAM use.
connDataFolder=fullfile(parentFolder, ...
    ['savedata_bivariate_' char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss-SSS'))]);
mkdir(connDataFolder);

addpath(fileparts(mfilename('fullpath')));
ft_defaults;
pool=gcp('nocreate');
if isempty(pool), pool=parpool('Processes',numWorkers); end
% Reuse an existing process pool; copy the active project/FieldTrip paths to workers.
f=parfevalOnAll(pool,@path,0,path); fetchOutputs(f);
f=parfevalOnAll(pool,@ft_defaults,0); fetchOutputs(f);

parfor (i=1:numel(subjects),numWorkers)
    saveConnDataGoodChannels_bivariate(subjects{i},protocolNameList, ...
        badEyeCondition,badTrialVersion,ftDataFolder,'granger',stRange,connDataFolder);
end
fprintf('Bivariate GC saved in: %s\n',connDataFolder);
