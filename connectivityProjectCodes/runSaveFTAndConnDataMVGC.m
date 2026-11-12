% MVGC1 copy of runSaveFTAndConnDataParametricFast.
% Fixed-order GC, with the same FieldTrip inputs and saved Pre/Post matrix layout.
% Install MVGC1 at add_params.mvgc_root; each worker initialises it on first use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connMethod = 'granger';
add_params = struct;
add_params.granger_model = 'bivariate';  % Fit each pair; 'multivariate' fits all selected channels.
add_params.mvgc_conditional = 'no';      % Joint model: 'yes' conditions on all other fitted channels.
% Joint + 'no' preserves nonconditional pair GC using autocovariance submatrices.
add_params.mvar_order = 20;             % Samples: 100 ms at fs=200 Hz.
add_params.resamplefs = 200;            % [] = original rate; downsampling only.
add_params.mvgc_regmode = 'OLS';        % MVGC regression: 'OLS' (QR) or 'LWR' (Morf).
add_params.mvgc_spectral_method = 'statespace'; % 'statespace' or 'autocov'; joint marginal uses autocov.
add_params.mvgc_demean = 'ensemble';     % MVGC pooled centering; 'pertrial' adds FieldTrip-like trial demeaning.
add_params.mvar_zscore = 'no';           % 'yes' = scale each channel across all trials.
add_params.mvar_freqstep = 1;            % Desired maximum frequency spacing in Hz.
add_params.mvgc_fres = [];               % [] = derive from freqstep; else number of intervals from 0 to Nyquist.
add_params.mvgc_acmaxlags = 0;           % Autocov path: 0 = automatic decay-based lag count.
add_params.mvgc_acdectol = 1e-8;         % Autocov path: decay tolerance.
add_params.remove_bad_eletrod_flag = true;
add_params.line_noise_freq = [50 100];   % Only lines strictly below the resampled Nyquist are removed.
add_params.ln_remove_method = 'zero';
add_params.mvgc_root = '/home/venkatesh/Documents/IISc/code/Meditation/MVGC1';
keepTrialsList = {'no'};
maxFrequencyList = 100;                 % Save up to this Hz; MVGC evaluates the full Nyquist range first.
lineNoiseRemovalList = true;
numWorkers = 8;                         % Subjects run independently; 0 = serial.

% Subjects, protocols, time windows and data locations.
useTheseIndices = [];                   % [] = all good subjects; e.g. 1:8 = the first eight good subjects.
protocolNameList = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2'};
stRange = [0.25 1.25];                   % Post window in seconds; Pre uses the same duration ending at zero.
badEyeCondition = 'ep';                 % 'ep' = eye-position rejection; 'wo' = without eye-position rejection.
badTrialVersion = 'v8';                 % Bad-trial version used in the input filenames.
folderSourceString = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/ftData';
parentFolder = '/home/venkatesh/Documents/IISc/data'; % Parent folder for the saved FC results.
saveFTDataFlag = 0;                     % 1 = create FieldTrip data first; 0 = use existing FieldTrip files.
saveConnDataFlag = 1;                   % 1 = compute and save FC; 0 = skip FC computation.

%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare the tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Capture the original inputs before the tests change any options.
inputParams = struct;
inputParams.connMethod = connMethod;
inputParams.add_params = add_params;
inputParams.keepTrialsList = keepTrialsList;
inputParams.maxFrequencyList = maxFrequencyList;
inputParams.lineNoiseRemovalList = lineNoiseRemovalList;
inputParams.numWorkers = numWorkers;
inputParams.useTheseIndices = useTheseIndices;
inputParams.protocolNameList = protocolNameList;
inputParams.stRange = stRange;
inputParams.badEyeCondition = badEyeCondition;
inputParams.badTrialVersion = badTrialVersion;
inputParams.folderSourceString = folderSourceString;
inputParams.parentFolder = parentFolder;
inputParams.saveFTDataFlag = saveFTDataFlag;
inputParams.saveConnDataFlag = saveConnDataFlag;


if strcmp(add_params.granger_model, 'bivariate')
    methodToRun = 'granger_mvgc_bivariate';
elseif strcmp(add_params.mvgc_conditional, 'yes')
    methodToRun = 'granger_mvgc_conditional';
else
    methodToRun = 'granger_mvgc_joint_marginal';
end
baseLineNoiseFreq = add_params.line_noise_freq;

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataRoot = fullfile(parentFolder, ['savedata_mvgc_' add_params.granger_model '_' add_params.mvgc_regmode '_' timestamp]);
ftDataFolder = folderSourceString;       % Location of the existing subject/protocol FieldTrip files.
goodSubjectList = getGoodSubjectsBK1;    % Subjects included in the analysis.
if isempty(useTheseIndices)
    useTheseIndices = 1:numel(goodSubjectList);
end
selectedSubjects = goodSubjectList(useTheseIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%% Save FieldTrip data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFTDataFlag
    [allSubjectNames, expDateList] = getDemographicDetails('BK1');
    for i = 1:numel(selectedSubjects)
        subjectName = selectedSubjects{i};
        disp(['Analyzing for subject ' subjectName]);
        expDate = expDateList{strcmp(subjectName, allSubjectNames)};
        saveFTData(subjectName, expDate, protocolNameList, folderSourceString, ...
            badEyeCondition, badTrialVersion, ftDataFolder);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Save connectivity data %%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveConnDataFlag
    % Trial-wise FC is deliberately left as a placeholder in this version.
    if any(strcmp(keepTrialsList, 'yes'))
        error('keeptrials = yes is a placeholder. Use keepTrialsList = {''no''} for now.');
    end

    % Run serially if zero workers are requested or the parallel toolbox is absent.
    if numWorkers > 0 && exist('parpool', 'file') ~= 0 ...
            && license('test', 'Distrib_Computing_Toolbox')
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('Processes', numWorkers);
        else
            numWorkers = min(numWorkers, pool.NumWorkers);
        end
    else
        numWorkers = 0;
    end

    for iFreq = 1:numel(maxFrequencyList)
        add_params.foilim_max = maxFrequencyList(iFreq);
        add_params.line_noise_freq = baseLineNoiseFreq;
        if ~isempty(add_params.resamplefs)
            add_params.line_noise_freq = baseLineNoiseFreq(baseLineNoiseFreq < add_params.resamplefs/2);
        end
        for iNoise = 1:numel(lineNoiseRemovalList)
            add_params.linenoise_flag = lineNoiseRemovalList(iNoise) && ~isempty(add_params.line_noise_freq);
            % Derive the label from the logical flag; a dictionary is not needed.
            if add_params.linenoise_flag
                lineNoiseLabel = 'withoutLineNoise';
            elseif lineNoiseRemovalList(iNoise)
                lineNoiseLabel = 'noLineFrequenciesInRange';
            else
                lineNoiseLabel = 'withLineNoise';
            end
            for iKeep = 1:numel(keepTrialsList)
                add_params.keeptrials = keepTrialsList{iKeep};
                testLabel = sprintf('p%d_%gHz_%s_keeptrials_%s', ...
                    add_params.mvar_order, add_params.foilim_max, lineNoiseLabel, add_params.keeptrials);
                connDataFolder = fullfile(connDataRoot, testLabel);
                disp(['Running test: ' methodToRun ' | ' testLabel]);

                % One shared parameter file per test, above the subject folders.
                testParams = struct;
                testParams.timestamp = timestamp;
                testParams.script = mfilename('fullpath');
                testParams.testLabel = testLabel;
                testParams.methodToRun = methodToRun;
                testParams.selectedSubjects = selectedSubjects;
                testParams.numWorkersUsed = numWorkers;
                testParams.ftDataFolder = ftDataFolder;
                testParams.connDataRoot = connDataRoot;
                testParams.connDataFolder = connDataFolder;
                testParams.add_params = add_params;
                saveTestInputParameters(connDataFolder, inputParams, testParams);

                % Each subject runs independently; zero workers means serial execution.
                parfor (i = 1:numel(selectedSubjects), numWorkers)
                    subjectName = selectedSubjects{i};
                    disp(['Analyzing for subject ' subjectName]);
                    pre_datetime = datetime;
                    saveConnDataWithGrangerGoodChannelsMVGC( ...
                        subjectName, protocolNameList, badEyeCondition, badTrialVersion, ...
                        ftDataFolder, methodToRun, stRange, connDataFolder, add_params);
                    dt = datetime - pre_datetime;
                    fprintf('Time taken: %d minutes and %.2f seconds\n', ...
                        floor(minutes(dt)), mod(seconds(dt), 60));
                end
            end
        end
    end
end

function saveTestInputParameters(connDataFolder, inputParams, testParams)
mkdir(connDataFolder);
fileID = fopen(fullfile(connDataFolder, 'input_parameters.txt'), 'w');
fileCleanup = onCleanup(@() fclose(fileID));
fprintf(fileID, 'MVGC connectivity: input parameters\n\n');
fprintf(fileID, 'Original sweep inputs\n%s\n\n', ...
    jsonencode(inputParams, 'PrettyPrint', true));
fprintf(fileID, 'Effective settings for this test\n%s\n', ...
    jsonencode(testParams, 'PrettyPrint', true));
end
