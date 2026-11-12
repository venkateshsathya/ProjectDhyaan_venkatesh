% Independent FC tests; subjects run in parallel within each test.
% Set the input options here. No options need to be edited in the second file.

%{
Created two separate files; both originals are unchanged:
- [runSaveFTAndConnDataMVGC.m](/home/venkatesh/Documents/IISc/code/Meditation/ProjectDhyaan_venkatesh/connectivityProjectCodes/runSaveFTAndConnDataMVGC.m)
- [saveConnDataWithGrangerGoodChannelsMVGC.m](/home/venkatesh/Documents/IISc/code/Meditation/ProjectDhyaan_venkatesh/connectivityProjectCodes/saveConnDataWithGrangerGoodChannelsMVGC.m)
Defaults are fs=200 Hz, p=20, pairwise nonconditional GC, OLS fitting, with your existing subjects, protocols, and eight-worker setting.
Install MVGC1 from an Ubuntu terminal:
git clone --depth 1 https://github.com/lcbarnett/MVGC1.git \
  /home/venkatesh/Documents/IISc/code/Meditation/MVGC1
Then run this from MATLAB with your connectivity code folder on its path:
runSaveFTAndConnDataMVGC
The helper initializes MVGC automatically on each worker. Use its startup.m for initialization rather than adding the entire directory tree with genpath. Official startup instructions
MVGC fits a vector autoregressive model, then derives time-domain and frequency-domain Granger causality from that model. Its fitting choices are OLS, using QR-based least squares, and LWR, using Morf's recursive algorithm. These are different implementation choices from your BioSig algorithm numbers. MVGC fitting documentation
Option	Your current FieldTrip code	New MVGC copy
VAR fitting	BioSig 1/2 or BSMART	mvgc_regmode = 'OLS' or 'LWR'
Model order	Fixed mvar_order	Same; defaults to 20
Spectral calculation	FieldTrip's VAR-to-spectrum pipeline	State-space default; autocovariance option
Centering	Per-trial demeaning	Pooled centering default; optional 'pertrial'
Joint GC	Currently explicitly nonconditional	Conditional or nonconditional
Fourier settings	DPSS/taper settings for the Fourier route	No Fourier multitaper estimation in this adapter


The state-space route is recommended in MVGC1's own demonstration as generally faster and potentially more accurate than its autocovariance route. Official demonstration
For preprocessing closer to your existing FieldTrip runs, set:
add_params.mvgc_demean = 'pertrial';
For joint conditional GC, controlling for the other fitted channels:
add_params.granger_model = 'multivariate';
add_params.mvgc_conditional = 'yes';
Leaving conditioning at 'no' preserves a nonconditional measure. Changing this option changes the scientific quantity being compared, not just the numerical algorithm.
MVGC's appeal for GC research is its focused workflow: model-order selection, model diagnostics, conditional spectral GC, and statistical inference. I haven't verified that it is more popular overall than FieldTrip. MVGC supports AIC/BIC selection and significance testing, but these copies retain fixed p and do not calculate p-values. Order-selection documentation, inference example
The copies save your usual Pre/Post matrices, shared parameter text files, time-domain GC, and model diagnostics. Synthetic tests passed for both regressions, conditioning, direction, resampling, and bad-channel masking. The full subject analysis has not been started.
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connMethod = 'granger';                  % 'granger', 'coh', 'ppc', 'plv', 'wpli', 'wpli_debiased'.
add_params.granger_method = 'parametric'; % 'parametric' or 'nonparametric'; used only for Granger.
add_params.granger_model = 'bivariate';   % 'bivariate' = fit each pair; 'multivariate' = fit all good channels.
% Multivariate refers to the fitted model. The GC measure remains nonconditional, as before.
keepTrialsList = {'no'};                 % 'no' = pooled FC; 'yes' is a placeholder and is not run yet.
maxFrequencyList = [100];                % Maximum analyzed frequency in Hz; e.g. [200 500]. Sampling rate is unchanged.
lineNoiseRemovalList = [true];           % true = removed; false = retained; [false true] runs both tests.
add_params.remove_bad_eletrod_flag = true; % true = exclude bad channels from calculation; false = include them.
% Bad rows/columns in the saved matrices remain NaN in both cases, as before.
numWorkers = 8;                         % 0 = serial; e.g. 8 = eight workers if the parallel toolbox is available.
add_params.resamplefs = 200; % [] = original sampling rate; 500 = resample to 500 Hz.
% Parametric Granger options; ignored by the nonparametric and phase methods.
add_params.mvar_order = 20;              % AR order in samples; 20 samples correspond to 20 ms at fs = 1000 Hz.
add_params.mvar_method = 'biosig';       % 'bsmart' or 'biosig'.
add_params.mvar_algorithm = 1;           % BioSig only: 1 = Yule-Walker; 2 = Vieira-Morf.
add_params.mvar_demean = 'yes';          % 'yes' or 'no'; subtract the temporal mean, not the average reference.
add_params.mvar_zscore = 'no';           % 'yes' or 'no'; standardize channel amplitudes before fitting.
add_params.mvar_taper = 'rectwin';       % 'rectwin' retains the existing AR fit without tapering.
add_params.mvar_freqstep = 1;            % Frequency grid spacing in Hz for evaluating the fitted AR model.

% Fourier options; used by nonparametric Granger, coherence, PPC, etc.
add_params.taper = 'dpss';               % 'dpss', 'hanning', 'hamming', or 'rectwin'.
add_params.tapsmofrq = 1;                % DPSS half-bandwidth in Hz; used only when taper = 'dpss'.

% Line-noise options; the same cleaning is available for every FC method.
add_params.line_noise_freq = [50 100];   % Frequencies to remove in Hz; e.g. [50 100 150].
add_params.ln_remove_method = 'zero';    % 'zero', 'neighbour', or 'neighbour_fft'.

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


% Use explicit method names in the saved filenames.
methodToRun = connMethod;
switch connMethod
    case 'granger'
        switch add_params.granger_method
            case 'parametric'
                switch add_params.granger_model
                    case 'bivariate'
                        methodToRun = 'granger_parametric';
                    case 'multivariate'
                        methodToRun = 'granger_parametric_multivariate';
                end
            case 'nonparametric'
                switch add_params.granger_model
                    case 'bivariate'
                        methodToRun = 'granger_nonparametric_bivariate';
                    case 'multivariate'
                        methodToRun = 'granger_nonparametric_multivariate';
                end
        end
end

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataRoot = fullfile(parentFolder, ['savedata_parametric_fast_' timestamp]);
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
        delete(gcp('nocreate'));
        parpool('Processes', numWorkers);
    else
        numWorkers = 0;
    end

    for iFreq = 1:numel(maxFrequencyList)
        add_params.foilim_max = maxFrequencyList(iFreq);
        for iNoise = 1:numel(lineNoiseRemovalList)
            add_params.linenoise_flag = lineNoiseRemovalList(iNoise);
            % Derive the label from the logical flag; a dictionary is not needed.
            if add_params.linenoise_flag
                lineNoiseLabel = 'withoutLineNoise';
            else
                lineNoiseLabel = 'withLineNoise';
            end
            for iKeep = 1:numel(keepTrialsList)
                add_params.keeptrials = keepTrialsList{iKeep};
                testLabel = sprintf('%dHz_%s_keeptrials_%s', ...
                    add_params.foilim_max, lineNoiseLabel, add_params.keeptrials);
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
                    saveConnDataWithGrangerGoodChannelsDFTParametricFast( ...
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
fprintf(fileID, 'Parametric fast connectivity: input parameters\n\n');
fprintf(fileID, 'Original sweep inputs\n%s\n\n', ...
    jsonencode(inputParams, 'PrettyPrint', true));
fprintf(fileID, 'Effective settings for this test\n%s\n', ...
    jsonencode(testParams, 'PrettyPrint', true));
end
