% Independent Granger tests; run subjects in parallel within each test.

maxFrequencyList = [200];
lineNoiseRemovalList = [true]; % false = with line noise; true = removed.
lineNoiseDict = dictionary(["true", "false"], ['withoutLineNoise', 'withLineNoise']);
lineNoiseLabels = repmat({''}, 1, numel(lineNoiseRemovalList));
idx = 1;
for lineNoiseval = lineNoiseRemovalList
    lineNoiseLabels{idx} = lineNoiseDict(lineNoiseval);
    idx = idx+1;
end
keepTrialsList = {'no'}; % Use {'no','yes'} to run pooled versus per-trial tests.
% numWorkers can be set to number of cores, I have used maximum of 8 on Ubuntu I9
% processor computer.
numWorkers =0; % Requires Parallel Computing Toolbox; use 0 for serial execution.
folderSourceString = '/run/media/venkatesh/Data1/data/FCDataFromVinodJiani/ftData';

% Location where FC data are stored.
parentFolder = '/home/venkatesh/Documents/IISc/data';

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataRoot = fullfile(parentFolder, ['savedata_parametric_fast_' timestamp]);

goodSubjectList = getGoodSubjectsBK1; % subjects for which data is computed
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

ftDataFolder = fullfile(folderSourceString); % fieldtrip data is saved here

badEyeCondition = 'ep'; % use 'wo' when eye data is not used and 'ep' for eye position. Note that bad trials will be removed before saving in fieldtrip format
badTrialVersion = 'v8';


%%%%%%%%%%%%%%%%%%%%%%%% Get data for all protocols %%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];

useTheseIndices = 1:length(goodSubjectList);
selectedSubjects = goodSubjectList(useTheseIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%% Save Fieldtrip data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFTDataFlag=0;
if saveFTDataFlag
    for i=1:length(useTheseIndices)
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
add_params.line_noise_freq = [50,100];
add_params.ln_remove_method = 'zero'; % Other methods are 'neighbour' or 'neighbour_fft'
stRange = [0.25 1.25];
saveConnDataFlag=1;
if saveConnDataFlag
    if numWorkers > 0 && exist('parpool','file') ~= 0 ...
            && license('test','Distrib_Computing_Toolbox')
        delete(gcp('nocreate'));
        parpool('Processes', numWorkers);
    else
        numWorkers = 0;   % Run serially.
    end
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
