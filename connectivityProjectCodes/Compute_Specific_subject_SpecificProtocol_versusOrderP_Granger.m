% Independent FC tests; subjects run in parallel within each test.
% Set the input options here. No options need to be edited in the second file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connMethod = 'granger';                  % 'granger', 'coh', 'ppc', 'plv', 'wpli', 'wpli_debiased'.
add_params.granger_method = 'parametric'; % 'parametric' or 'nonparametric'; used only for Granger.
add_params.granger_model = 'bivariate';   % 'bivariate' = fit each pair; 'multivariate' = fit all good channels.
% Multivariate refers to the fitted model. The GC measure remains nonconditional, as before.
keepTrialsList = {'no'};                 % 'no' = pooled FC; 'yes' is a placeholder and is not run yet.
maxFrequencyList = [200];                % Maximum analyzed frequency in Hz; e.g. [200 500]. Sampling rate is unchanged.
lineNoiseRemovalList = [true];           % true = removed; false = retained; [false true] runs both tests.
add_params.remove_bad_eletrod_flag = true; % true = exclude bad channels from calculation; false = include them.
% Bad rows/columns in the saved matrices remain NaN in both cases, as before.
numWorkers = 0;                         % 0 = serial; e.g. 8 = eight workers if the parallel toolbox is available.
add_params.resamplefs = [200];%200; % [] = original sampling rate; 500 = resample to 500 Hz.

% Parametric Granger options; ignored by the nonparametric and phase methods.
mvarOrderList = [5 10 20 40 80 160 320];            % AR orders p to compare for parametric Granger.
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
protocolNameList = {'G1','G2', 'M1', 'M2'};            % Test one protocol; change to {'G1'}, {'M1'}, etc.
stRange = [0.25 1.25];                   % Post window in seconds; Pre uses the same duration ending at zero.
badEyeCondition = 'ep';                 % 'ep' = eye-position rejection; 'wo' = without eye-position rejection.
badTrialVersion = 'v8';                 % Bad-trial version used in the input filenames.
folderSourceString = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/ftData';
parentFolder = '/home/venkatesh/Documents/IISc/data'; % Parent folder for the saved FC results.
saveFTDataFlag = 0;                     % 1 = create FieldTrip data first; 0 = use existing FieldTrip files.
saveConnDataFlag = 1;                   % 1 = compute and save FC; 0 = skip FC computation.

%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare the tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep the existing method names in saved filenames, for downstream plotting.
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
                        methodToRun = 'granger_bivariate';
                    case 'multivariate'
                        methodToRun = 'granger';
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
if ~isempty(selectedSubjectsOverride)
    selectedSubjects = selectedSubjectsOverride;
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

    for iOrder = 1:numel(mvarOrderList)             % Repeat the same test for each AR order p.
        add_params.mvar_order = mvarOrderList(iOrder);
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
                    testLabel = sprintf('p%d_%dHz_%s_keeptrials_%s', ...
                        add_params.mvar_order, add_params.foilim_max, lineNoiseLabel, add_params.keeptrials);
                    connDataFolder = fullfile(connDataRoot, testLabel);
                    disp(['Running test: ' methodToRun ' | ' testLabel]);

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
        end % iFreq
    end % iOrder
end