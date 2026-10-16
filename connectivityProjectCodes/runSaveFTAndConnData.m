% Save data in fieldtrip format

% Mandatory fixed options
%folderSourceString = '/mnt/NeoLabData/Projects/ProjectDhyaan/BK1/';%'N:\Projects\ProjectDhyaan\BK1';
folderSourceString = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1';

% Location where FC data are stored.
parentFolder = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1';

timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
connDataFolder = fullfile(parentFolder, ['savedata_' timestamp]);
% folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\ProjectDhyaan\BK1'; % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData

goodSubjectList = getGoodSubjectsBK1; % subjects for which data is computed
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

ftDataFolder = fullfile(folderSourceString,'data','ftData'); % fieldtrip data is saved here
%makeDirectory(ftDataFolder);

badEyeCondition = 'ep'; % use 'wo' when eye data is not used and 'ep' for eye position. Note that bad trials will be removed before saving in fieldtrip format
badTrialVersion = 'v8';


%%%%%%%%%%%%%%%%%%%%%%%% Get data for all protocols %%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 

useTheseIndices = 1:length(goodSubjectList);

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
% connMethod options: granger_bivariate, 'granger_multivariate','granger'
connMethod = 'granger_bivariate';
% global remove_bad_eletrod_flag
add_params.remove_bad_eletrod_flag = true;
% global foilim_max
add_params.foilim_max = 200;
% global linenoise_flag
add_params.linenoise_flag = true;
% global line_noise_freq 
add_params.line_noise_freq = [50,100];
% global ln_remove_method
add_params.ln_remove_method = 'zero'; % Other methods are 'neighbour' or 'neighbour_fft'
% fla
stRange = [0.25 1.25];
saveConnDataFlag=1;
if saveConnDataFlag
    for i=1:length(useTheseIndices) %#ok<*UNRCH>
        subjectName = goodSubjectList{useTheseIndices(i)};
        disp(['Analyzing for subject ' subjectName]);
        expDate = expDateList{strcmp(subjectName,allSubjectNames)};
        pre_datetime = datetime;
        saveConnDataWithGrangerGoodChannelsDFT(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange, connDataFolder, add_params);
        % fprintf("FC ", connMethod, " bivariate")
        post_datetime = datetime;
        dt = post_datetime - pre_datetime;
        fprintf('Time taken: %d minutes and %.2f seconds\n', ...
            floor(minutes(dt)), mod(seconds(dt),60));
    end
end 
