% Independent Granger tests; run subjects in parallel within each test.
% 200/500 Hz are maximum analyzed frequencies; sampling rate is unchanged.
maxFrequencyList = [200];
lineNoiseRemovalList = [false true]; % false = with line noise; true = removed.
lineNoiseLabels = {'withoutLineNoise'};
keepTrialsList = {'yes'}; % Use {'no','yes'} to run pooled and per-trial tests.
numWorkers = 8; % Requires Parallel Computing Toolbox; use 0 for serial execution.

% Mandatory fixed options
%folderSourceString = '/mnt/NeoLabData/Projects/ProjectDhyaan/BK1/';%'N:\Projects\ProjectDhyaan\BK1';
folderSourceString = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/ftData';

% Location where FC data are stored.
parentFolder = '/home/venkatesh/Documents/IISc/data';

existingFCFolder = '';%'/home/venkatesh/Documents/IISc/data/savedata_parametric_fast_2026-10-23_20-55-34/500Hz_withLineNoise_keeptrials_no'; % Optional: previous run root OR one test folder; '' starts a new run.
% Keep the AR settings, time windows and DFT removal method matched to the original run.

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataRoot = fullfile(parentFolder, ['savedata_parametric_fast_' timestamp]);
if ~isempty(existingFCFolder)
    existingFCFolder = regexprep(existingFCFolder, '[\\/]+$', ''); % Remove trailing separator.
    [runRoot,testFolder] = fileparts(existingFCFolder);
    if endsWith(testFolder, {'_keeptrials_no','_keeptrials_yes'})
        % One test folder: resume only its frequency, line-noise and trial settings.
        connDataRoot = runRoot;
        testParts = strsplit(testFolder, '_');
        maxFrequencyList = str2double(erase(testParts{1}, 'Hz'));
        lineNoiseLabels = testParts(2);
        lineNoiseRemovalList = strcmp(testParts{2}, 'withoutLineNoise');
        keepTrialsList = testParts(4);
    else
        connDataRoot = existingFCFolder; % Resume all combinations selected above.
    end
end
% folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\ProjectDhyaan\BK1'; % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData

goodSubjectList = getGoodSubjectsBK1; % subjects for which data is computed
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

ftDataFolder = fullfile(folderSourceString); % fieldtrip data is saved here
%makeDirectory(ftDataFolder);

badEyeCondition = 'ep'; % use 'wo' when eye data is not used and 'ep' for eye position. Note that bad trials will be removed before saving in fieldtrip format
badTrialVersion = 'v8';


%%%%%%%%%%%%%%%%%%%%%%%% Get data for all protocols %%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];

useTheseIndices = 1:length(goodSubjectList);
selectedSubjects = goodSubjectList(useTheseIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%% Save Fieldtrip data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFTDataFlag=0;
if saveFTDataFlag
    for i=1:length(useTheseIndices) %#ok<*UNRCH>
        subjectName = goodSubjectList{useTheseIndices(i)};
        disp(['Analyzing for subject ' subjectName]);
        expDate = expDateList{strcmp(subjectName,allSubjectNames)};

        saveFTData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,ftDataFolder);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% Save Connectivity data %%%%%%%%%%%%%%%%%%%%%%%%%%%
% connMethod options: 'granger_parametric', 'granger_bivariate', 'granger_multivariate', 'granger'
connMethod = 'granger_parametric';
% AR settings: order is illustrative and should be selected/validated for
% the data.MVAR: multivariate autogressive
add_params.mvar_order = 20;        % Samples: 20 ms at fs = 1000 Hz
add_params.mvar_method = 'bsmart'; % 'bsmart' or 'biosig'
add_params.mvar_algorithm = 1;     % BioSig only: 2 = Vieira-Morf, 1 = Yule-Walker
add_params.mvar_demean = 'yes';    % Temporal mean removal; not average referencing
% global remove_bad_eletrod_flag
add_params.remove_bad_eletrod_flag = true;
% global foilim_max
% add_params.foilim_max is set separately for each test below.
% global linenoise_flag
% add_params.linenoise_flag is set separately for each test below.
% global line_noise_freq
add_params.line_noise_freq = [50,100];
% global ln_remove_method
add_params.ln_remove_method = 'zero'; % Other methods are 'neighbour' or 'neighbour_fft'
% fla
stRange = [0.25 1.25];
saveConnDataFlag=1;
if saveConnDataFlag
    delete(gcp('nocreate'));
    parpool('Processes', numWorkers);
    for iFreq = 1:numel(maxFrequencyList)
        add_params.foilim_max = maxFrequencyList(iFreq);
        for iNoise = 1:numel(lineNoiseRemovalList)
            add_params.linenoise_flag = lineNoiseRemovalList(iNoise);
            for iKeep = 1:numel(keepTrialsList)
                add_params.keeptrials = keepTrialsList{iKeep};
                testLabel = sprintf('%dHz_%s_keeptrials_%s', ...
                    add_params.foilim_max, lineNoiseLabels{iNoise}, add_params.keeptrials);
                connDataFolder = fullfile(connDataRoot, testLabel);
                disp(['Running test: ' testLabel]);
                parfor (i = 1:numel(selectedSubjects), numWorkers)
                    subjectName = selectedSubjects{i};
                    disp(['Analyzing for subject ' subjectName]);
                    pre_datetime = datetime;
                    saveConnDataWithGrangerGoodChannelsDFTParametricFast(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange, connDataFolder, add_params);
                    % fprintf("FC ", connMethod, " bivariate")
                    post_datetime = datetime;
                    dt = post_datetime - pre_datetime;
                    fprintf('Time taken: %d minutes and %.2f seconds\n', ...
                        floor(minutes(dt)), mod(seconds(dt),60));
                end
            end
        end
    end
end
