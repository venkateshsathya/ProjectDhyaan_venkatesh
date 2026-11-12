% Independent sample-rate x AR-order FC tests; subjects run in parallel.
% Copy of Compute_Specific_subject_SpecificProtocol_versusOrderP_Granger.
% Run with the original connectivityProjectCodes folder and FieldTrip on the MATLAB path.
% Every rate starts from the original FieldTrip input, never a previous resample.
% Set the input options here. No options need to be edited in the second file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connMethod = 'granger';                  % 'granger', 'coh', 'ppc', 'plv', 'wpli', 'wpli_debiased'.
add_params = struct;                    % Reset options when rerunning the script.
add_params.granger_method = 'parametric'; % 'parametric' or 'nonparametric'; used only for Granger.
add_params.granger_model = 'bivariate';   % 'bivariate' = fit each pair; 'multivariate' = fit all good channels.
% Multivariate refers to the fitted model. The GC measure remains nonconditional, as before.
keepTrialsList = {'no'};                 % 'no' = pooled FC; 'yes' is a placeholder and is not run yet.
lineNoiseRemovalList = false;             % true = removed; false = retained; [false true] runs both tests.
add_params.remove_bad_eletrod_flag = true; % true = exclude bad channels from calculation; false = include them.
% Bad rows/columns in the saved matrices remain NaN in both cases, as before.
numWorkers = 0;                         % 0 = serial; e.g. 8 = eight workers if the parallel toolbox is available.
sampleRateList = [200];%[20 50 100 200];      % Target rates in Hz; use rates no higher than the input rate.
maxFrequencyList = sampleRateList / 2;   % One Nyquist limit per corresponding sample rate.

% Parametric Granger options; ignored by the nonparametric and phase methods.
mvarOrderList = [49 99 199  249];            % AR orders p to compare for parametric Granger.
add_params.mvar_order = 20;              % Default AR order; overwritten by mvarOrderList during the sweep.
add_params.mvar_method = 'bsmart';       % 'bsmart' or 'biosig'.
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
selectedSubjectsOverride = {'053DR'};    % {'053DR'} = test one subject; [] = use useTheseIndices.
protocolNameList = {'G1','G2', 'M1', 'M2'}; % Protocols to test; use {'G1'} for one protocol.
stRange = [0.25 1.25];                   % Post window in seconds; Pre uses the same duration ending at zero.
badEyeCondition = 'ep';                 % 'ep' = eye-position rejection; 'wo' = without eye-position rejection.
badTrialVersion = 'v8';                 % Bad-trial version used in the input filenames.
folderSourceString = '/home/venkatesh/Documents/IISc/data/Vinod_SegmentedFTData/ftData';
parentFolder = '/home/venkatesh/Documents/IISc/data'; % Parent folder for the saved FC results.
saveFTDataFlag = 0;                     % 1 = create FieldTrip data first; 0 = use existing FieldTrip files.
saveConnDataFlag = 1;                   % 1 = compute and save FC; 0 = skip FC computation.

%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare the tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Capture the original inputs before the sweep changes any options.
inputParams = struct;
inputParams.connMethod = connMethod;
inputParams.add_params = add_params;
inputParams.keepTrialsList = keepTrialsList;
inputParams.lineNoiseRemovalList = lineNoiseRemovalList;
inputParams.numWorkers = numWorkers;
inputParams.sampleRateList = sampleRateList;
inputParams.maxFrequencyList = maxFrequencyList;
inputParams.mvarOrderList = mvarOrderList;
inputParams.useTheseIndices = useTheseIndices;
inputParams.selectedSubjectsOverride = selectedSubjectsOverride;
inputParams.protocolNameList = protocolNameList;
inputParams.stRange = stRange;
inputParams.badEyeCondition = badEyeCondition;
inputParams.badTrialVersion = badTrialVersion;
inputParams.folderSourceString = folderSourceString;
inputParams.parentFolder = parentFolder;
inputParams.saveFTDataFlag = saveFTDataFlag;
inputParams.saveConnDataFlag = saveConnDataFlag;

baseLineNoiseFreq = add_params.line_noise_freq;

% Use the method names recognized by the calculation helper.
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
isParametricGranger = strcmp(connMethod, 'granger') ...
    && strcmp(add_params.granger_method, 'parametric');
ordersToRun = mvarOrderList;
if ~isParametricGranger
    ordersToRun = NaN;                  % These methods have no AR order; run once per rate.
end

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataRoot = fullfile(parentFolder, ['savedata_granger_fs_order_' timestamp]);
ftDataFolder = folderSourceString;       % Location of the existing subject/protocol FieldTrip files.
if ~isempty(selectedSubjectsOverride)
    selectedSubjects = selectedSubjectsOverride;
else
    goodSubjectList = getGoodSubjectsBK1;
    if isempty(useTheseIndices)
        useTheseIndices = 1:numel(goodSubjectList);
    end
    selectedSubjects = goodSubjectList(useTheseIndices);
end

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

    for iRate = 1:numel(sampleRateList)
        add_params.resamplefs = sampleRateList(iRate);
        add_params.foilim_max = maxFrequencyList(iRate);
        % Remove only line frequencies strictly below this rate's Nyquist limit.
        add_params.line_noise_freq = baseLineNoiseFreq(baseLineNoiseFreq < add_params.foilim_max);

        ordersThisRate = ordersToRun;
        if isParametricGranger
            % Fit only orders below the number of samples in each time window.
            ordersThisRate = ordersToRun(ordersToRun < floor(diff(stRange) * add_params.resamplefs));
        end
        for iOrder = 1:numel(ordersThisRate)
            orderLabel = 'pNA';
            if isParametricGranger
                add_params.mvar_order = ordersThisRate(iOrder);
                orderLabel = sprintf('p%d', add_params.mvar_order);
            end
            for iNoise = 1:numel(lineNoiseRemovalList)
                % An empty list makes FieldTrip default to 50 Hz; disable the filter instead.
                add_params.linenoise_flag = lineNoiseRemovalList(iNoise) ...
                    && ~isempty(add_params.line_noise_freq);
                if add_params.linenoise_flag
                    lineNoiseLabel = 'withoutLineNoise';
                elseif lineNoiseRemovalList(iNoise)
                    lineNoiseLabel = 'noLineFrequenciesInRange';
                else
                    lineNoiseLabel = 'withLineNoise';
                end
                for iKeep = 1:numel(keepTrialsList)
                    add_params.keeptrials = keepTrialsList{iKeep};
                    testLabel = sprintf('fs%.17gHz_%s_fmax%.17gHz_%s_keeptrials_%s', ...
                        add_params.resamplefs, orderLabel, add_params.foilim_max, ...
                        lineNoiseLabel, add_params.keeptrials);
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
                    if numWorkers == 0
                        for i = 1:numel(selectedSubjects)
                            runSweepSubject(selectedSubjects{i}, protocolNameList, ...
                                badEyeCondition, badTrialVersion, ftDataFolder, ...
                                methodToRun, stRange, connDataFolder, add_params);
                        end
                    else
                        parfor (i = 1:numel(selectedSubjects), numWorkers)
                            runSweepSubject(selectedSubjects{i}, protocolNameList, ...
                                badEyeCondition, badTrialVersion, ftDataFolder, ...
                                methodToRun, stRange, connDataFolder, add_params);
                        end
                    end
                end
            end
        end % iOrder
    end % iRate
end

function runSweepSubject(subjectName, protocolNameList, badEyeCondition, ...
    badTrialVersion, ftDataFolder, methodToRun, stRange, connDataFolder, add_params)
disp(['Analyzing for subject ' subjectName]);
pre_datetime = datetime;
saveConnDataWithGrangerGoodChannelsDFTParametricFast( ...
    subjectName, protocolNameList, badEyeCondition, badTrialVersion, ...
    ftDataFolder, methodToRun, stRange, connDataFolder, add_params);
dt = datetime - pre_datetime;
fprintf('Time taken: %d minutes and %.2f seconds\n', ...
    floor(minutes(dt)), mod(seconds(dt), 60));
end

function saveTestInputParameters(connDataFolder, inputParams, testParams)
mkdir(connDataFolder);
fileID = fopen(fullfile(connDataFolder, 'input_parameters.txt'), 'w');
fileCleanup = onCleanup(@() fclose(fileID));
fprintf(fileID, 'Granger sample-rate/order sweep: input parameters\n\n');
fprintf(fileID, 'Original sweep inputs\n%s\n\n', ...
    jsonencode(inputParams, 'PrettyPrint', true));
fprintf(fileID, 'Effective settings for this test\n%s\n', ...
    jsonencode(testParams, 'PrettyPrint', true));
end
