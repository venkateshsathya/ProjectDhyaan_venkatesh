% Parallelise independent BioSig sampling-rate / AR-order combinations.
% Each worker owns one complete test setting and loops over subjects/protocols.
% There is no nested subject parfor. Both experiments share a queue of jobs.
% Uses the existing saveConnDataWithGrangerGoodChannelsDFTParametricFast helper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numWorkers = 8;                         % Concurrent processes; 0 = serial.
mvarAlgorithmList = 2;                  % Current BioSig option; e.g. [1 2 5 6 7] to compare several.

% Test 1: fixed p, varied sampling rate.
test1Order = 60;
test1SampleRates = [50 100 200 300 400 500 600 700 800 900 1000];
% Test 2: each requested order at both sampling rates.
test2SampleRates = [] ; % [200 500];
test2Orders =[];% [5 10 25 50 100 150 200 250 300 350 400 450 500];

add_params = struct;
add_params.granger_method = 'parametric';
add_params.granger_model = 'bivariate';  % Or 'multivariate'; GC remains nonconditional.
add_params.mvar_method = 'biosig';
add_params.mvar_demean = 'yes';
add_params.mvar_zscore = 'no';
add_params.mvar_taper = 'rectwin';
add_params.mvar_freqstep = 1;
add_params.keeptrials = 'no';
add_params.remove_bad_eletrod_flag = true;
add_params.line_noise_freq = [50 100];   % Applied only strictly below the resampled Nyquist limit.
add_params.ln_remove_method = 'zero';
add_params.taper = 'dpss';
add_params.tapsmofrq = 1;
lineNoiseRemoval = true;
% Each job sets resamplefs=Fs, mvar_order=p and foilim_max=Fs/2 automatically.

useTheseIndices = [];                   % [] = all good subjects.
selectedSubjectsOverride = {'053DR'};          % E.g. {'053DR'} for a one-subject smoke test.
protocolNameList = {'G2', 'M2'};%{'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
% protocolNameList = {'G1'};             % One protocol for a smoke test.
stRange = [0.25 1.25];                   % Post window; Pre has the same duration, ending at zero.
badEyeCondition = 'ep';
badTrialVersion = 'v8';
folderSourceString = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/ftData';
parentFolder = '/home/venkatesh/Documents/IISc/data';

%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare the jobs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputParams = struct;
inputParams.numWorkers = numWorkers;
inputParams.mvarAlgorithmList = mvarAlgorithmList;
inputParams.test1Order = test1Order;
inputParams.test1SampleRates = test1SampleRates;
inputParams.test2SampleRates = test2SampleRates;
inputParams.test2Orders = test2Orders;
inputParams.add_params = add_params;
inputParams.lineNoiseRemoval = lineNoiseRemoval;
inputParams.useTheseIndices = useTheseIndices;
inputParams.selectedSubjectsOverride = selectedSubjectsOverride;
inputParams.protocolNameList = protocolNameList;
inputParams.stRange = stRange;
inputParams.badEyeCondition = badEyeCondition;
inputParams.badTrialVersion = badTrialVersion;
inputParams.folderSourceString = folderSourceString;
inputParams.parentFolder = parentFolder;

if ~isempty(selectedSubjectsOverride)
    selectedSubjects = selectedSubjectsOverride;
else
    goodSubjectList = getGoodSubjectsBK1;
    if isempty(useTheseIndices), useTheseIndices = 1:numel(goodSubjectList); end
    selectedSubjects = goodSubjectList(useTheseIndices);
end
if strcmp(add_params.granger_model,'bivariate')
    methodToRun = 'granger_parametric';
else
    methodToRun = 'granger_parametric_multivariate';
end

timestamp = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss-SSS'));
connDataRoot = fullfile(parentFolder,['savedata_biosig_parallel_sweeps_' timestamp]);
mkdir(connDataRoot);

% Keep the requested settings in the manifest, including unfit orders.
jobs = struct('testName',{},'algorithm',{},'sampleRate',{},'order',{}, ...
    'folder',{},'skipReason',{});
for algorithm = mvarAlgorithmList(:)'
    for fs = test1SampleRates(:)'
        jobs(end+1) = makeJob('test1_fixed_order',algorithm,fs,test1Order,stRange, ...
            connDataRoot,add_params,lineNoiseRemoval); %#ok<SAGROW>
    end
    for fs = test2SampleRates(:)'
        for p = test2Orders(:)'
            jobs(end+1) = makeJob('test2_order_sweep',algorithm,fs,p,stRange, ...
                connDataRoot,add_params,lineNoiseRemoval); %#ok<SAGROW>
        end
    end
end
runnable = find(cellfun(@isempty,{jobs.skipReason}));
numWorkersUsed = 0;
if numWorkers > 0 && ~isempty(runnable) && exist('parpool','file') ~= 0 ...
        && license('test','Distrib_Computing_Toolbox')
    requestedWorkers = min(numWorkers,numel(runnable));
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers < requestedWorkers
        delete(pool); pool = [];
    end
    if isempty(pool), pool = parpool('Processes',requestedWorkers); end
    numWorkersUsed = min(requestedWorkers,pool.NumWorkers);
end

config = struct;
config.add_params = add_params;
config.lineNoiseRemoval = lineNoiseRemoval;
config.subjects = selectedSubjects;
config.protocolNameList = protocolNameList;
config.stRange = stRange;
config.badEyeCondition = badEyeCondition;
config.badTrialVersion = badTrialVersion;
config.ftDataFolder = folderSourceString;
config.methodToRun = methodToRun;
config.connDataRoot = connDataRoot;
config.timestamp = timestamp;
config.script = mfilename('fullpath');
config.numWorkersUsed = numWorkersUsed;

summary = repmat(emptySummary(),1,numel(jobs));
for iJob = 1:numel(jobs)
    summary(iJob) = describeJob(jobs(iJob));
    if ~isempty(jobs(iJob).skipReason)
        summary(iJob).status = 'skipped';
        summary(iJob).message = jobs(iJob).skipReason;
    end
end
% A complete manifest exists before any expensive connectivity calculation.
writetable(struct2table(summary,'AsArray',true),fullfile(connDataRoot,'job_summary.csv'));
writeParameters(fullfile(connDataRoot,'input_parameters.txt'),inputParams,config);
fprintf('Requested %d settings; %d runnable; %d skipped. Workers used: %d.\n', ...
    numel(jobs),numel(runnable),numel(jobs)-numel(runnable),numWorkersUsed);
fprintf('Results: %s\n',connDataRoot);

runnableJobs = jobs(runnable);
results = repmat(emptySummary(),1,numel(runnable));
if numWorkersUsed == 0
    for iRun = 1:numel(runnable)
        results(iRun) = runOneJob(runnableJobs(iRun),config,inputParams);
    end
else
    % A free worker takes the next (algorithm, Fs, p) setting automatically.
    parfor (iRun = 1:numel(runnable),numWorkersUsed)
        results(iRun) = runOneJob(runnableJobs(iRun),config,inputParams);
    end
end
summary(runnable) = results;
writetable(struct2table(summary,'AsArray',true),fullfile(connDataRoot,'job_summary.csv'));
fprintf('Finished: %d completed, %d with failures, %d skipped.\n', ...
    sum(strcmp({summary.status},'completed')), ...
    sum(strcmp({summary.status},'partial_failure') | strcmp({summary.status},'failed')), ...
    sum(strcmp({summary.status},'skipped')));

function job = makeJob(testName,algorithm,fs,p,stRange,root,options,removeLineNoise)
linesInRange = options.line_noise_freq(options.line_noise_freq < fs/2);
if removeLineNoise && ~isempty(linesInRange)
    noiseLabel = 'withoutLineNoise';
elseif removeLineNoise
    noiseLabel = 'noLineFrequenciesInRange';
else
    noiseLabel = 'withLineNoise';
end
label = sprintf('fs%gHz_p%d_fmax%gHz_%s_keeptrials_no',fs,p,fs/2,noiseLabel);
job = struct('testName',testName,'algorithm',algorithm,'sampleRate',fs,'order',p, ...
    'folder',fullfile(root,testName,sprintf('biosig%d',algorithm),label),'skipReason','');
if p >= floor(diff(stRange)*fs)
    job.skipReason = sprintf('p=%d >= nominal samples/window=%d (window=%g s, Fs=%g Hz).', ...
        p,floor(diff(stRange)*fs),diff(stRange),fs);
end
end

function result = runOneJob(job,config,inputParams)
result = describeJob(job);
options = config.add_params;
options.mvar_algorithm = job.algorithm;
options.resamplefs = job.sampleRate;
options.mvar_order = job.order;
options.foilim_max = job.sampleRate/2;
options.line_noise_freq = options.line_noise_freq(options.line_noise_freq < job.sampleRate/2);
options.linenoise_flag = config.lineNoiseRemoval && ~isempty(options.line_noise_freq);
mkdir(job.folder);
effective = config;
effective.job = job;
effective.add_params = options;
writeParameters(fullfile(job.folder,'input_parameters.txt'),inputParams,effective);
logID = fopen(fullfile(job.folder,'execution_log.txt'),'w');
closeLog = onCleanup(@() fclose(logID));
fprintf(logID,'Started %s | BioSig %d | Fs=%g Hz | p=%d\n', ...
    char(datetime('now')),job.algorithm,job.sampleRate,job.order);
fprintf('Starting %s | BioSig %d | Fs=%g Hz | p=%d\n', ...
    job.testName,job.algorithm,job.sampleRate,job.order);
startTime = tic;
for iSubject = 1:numel(config.subjects)
    subjectName = config.subjects{iSubject};
    for iProtocol = 1:numel(config.protocolNameList)
        protocolName = config.protocolNameList{iProtocol};
        try
            saveConnDataWithGrangerGoodChannelsDFTParametricFast( ...
                subjectName,{protocolName},config.badEyeCondition,config.badTrialVersion, ...
                config.ftDataFolder,config.methodToRun,config.stRange,job.folder,options);
            result.completedSubjectProtocols = result.completedSubjectProtocols + 1;
            fprintf(logID,'COMPLETED %s %s\n',subjectName,protocolName);
        catch runError
            % Preserve other subjects/protocols and other jobs when one calculation fails.
            result.failedSubjectProtocols = result.failedSubjectProtocols + 1;
            if isempty(result.message), result.message = runError.message; end
            fprintf(logID,'FAILED %s %s\n%s\n',subjectName,protocolName, ...
                getReport(runError,'extended','hyperlinks','off'));
        end
    end
end
result.elapsedSeconds = toc(startTime);
if result.failedSubjectProtocols == 0
    result.status = 'completed';
elseif result.completedSubjectProtocols == 0
    result.status = 'failed';
else
    result.status = 'partial_failure';
end
fprintf(logID,'Finished %s | %s | completed=%d failed=%d elapsed=%g s\n', ...
    char(datetime('now')),result.status,result.completedSubjectProtocols, ...
    result.failedSubjectProtocols,result.elapsedSeconds);
% A job-owned status file survives even if the client stops before final aggregation.
writetable(struct2table(result,'AsArray',true),fullfile(job.folder,'job_status.csv'));
fprintf('Finished %s | BioSig %d | Fs=%g Hz | p=%d | %s\n', ...
    job.testName,job.algorithm,job.sampleRate,job.order,result.status);
end

function result = emptySummary()
result = struct('testName','','algorithm',0,'sampleRate',0,'order',0, ...
    'status','pending','completedSubjectProtocols',0,'failedSubjectProtocols',0, ...
    'elapsedSeconds',0,'message','','outputFolder','');
end

function result = describeJob(job)
result = emptySummary();
result.testName = job.testName;
result.algorithm = job.algorithm;
result.sampleRate = job.sampleRate;
result.order = job.order;
result.outputFolder = job.folder;
end

function writeParameters(fileName,inputParams,effective)
fileID = fopen(fileName,'w');
closeFile = onCleanup(@() fclose(fileID));
fprintf(fileID,'BioSig parallel sweeps: input parameters\n\n');
fprintf(fileID,'Original inputs\n%s\n\n',jsonencode(inputParams,'PrettyPrint',true));
fprintf(fileID,'Effective settings\n%s\n',jsonencode(effective,'PrettyPrint',true));
end
