
function saveConnDataWithGranger(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange, connDataFolder)
% connMethod: 'granger_bivariate', 'granger_multivariate', or an existing
% FieldTrip method such as 'ppc' or 'coh'. 'granger' uses multivariate GC.

if ~exist('stRange','var');         stRange = [0.25 1.25];              end

numProtocols = length(protocolNameList);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder);
end
end

function saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder)

%connDataFolder = 'savedData'; % local project directory
%makeDirectory(connDataFolder);


[success, message] = mkdir(connDataFolder);
if ~success
    error('Could not create folder: %s', message);
end

disp(connDataFolder);
connDataFolder = fullfile(connDataFolder,subjectName);
makeDirectory(connDataFolder);

% Analysis file - this is where the final analyzed data is saved
analysisDetailsFileConn = fullfile(connDataFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']);

% ftData - data in fieldtrip format is saved in ftData
ftDataFileName = fullfile(ftDataFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
tmpData = load(ftDataFileName);
numGoodTrials = tmpData.numGoodTrials;

if numGoodTrials>0    
    data = tmpData.data;

    cfg        = [];
    cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0]; % just to have data of same length in Pre, Post


    dataForFT_NoTimevals = data;
    % if isfield(dataForFT_NoTimevals, 'timeVals')
    %     dataForFT_NoTimevals = rmfield(dataForFT_NoTimevals, 'timeVals');
    % end

    dataPre = ft_redefinetrial(cfg, dataForFT_NoTimevals);

    %dataPre   = ft_redefinetrial(cfg,data);
    
    cfg        = [];
    cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
    dataPost   = ft_redefinetrial(cfg,dataForFT_NoTimevals);
    
    [connPre,freqPre] = getConnIndividualSubject(dataPre,connMethod);
    [connPost,freqPost] = getConnIndividualSubject(dataPost,connMethod);
else
    connPre=[]; connPost=[]; freqPre=[]; freqPost=[];
end
save(analysisDetailsFileConn,'connPre','connPost','freqPre','freqPost','numGoodTrials');

end

function [conn,freqVals] = getConnIndividualSubject(data,connMethod)

[temp,fieldtripMethod]=getConnThisCondition(data,connMethod);
if strcmp(connMethod,'granger_bivariate')
    % Granger bivariate provides the output in a fashion of NC2*2 -
    % pairwise selected Granger FC and into 2 for two directions. 
    % We need to move data to a format of 64 by 64 by 2.
    % Extracts the pair of labels from temp.labelcmb, first column of pairLabels is
    % source, second is target.
    pairLabels = regexprep(temp.labelcmb,'\[.*$','');
    %Converts the labels in text such as Oz etc to numbers.
    [~,sourceIdx] = ismember(pairLabels(:,1),data.label);
    [~,targetIdx] = ismember(pairLabels(:,2),data.label);
    conn = zeros(numel(data.label),numel(data.label),numel(temp.freq));
    for k=1:size(pairLabels,1)
        conn(sourceIdx(k),targetIdx(k),:) = reshape(abs(temp.grangerspctrm(k,:)),1,1,[]);
    end
else
    conn = abs(temp.([fieldtripMethod 'spctrm']));
end
badElecs = data.badElecs;

conn(badElecs,:,:) = NaN;
conn(:,badElecs,:) = NaN;

freqVals = temp.freq;
end
function [conn_stat,fieldtripMethod] = getConnThisCondition(data,method)
cfg              = [];
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;
if(strcmp(method,'coh'))
    cfg.output = 'powandcsd';
else
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
end
cfg.foilim       = [0 200];
freqPost = ft_freqanalysis(cfg, data);

if(strcmp(method,'coh'))
    TfreqPost = ft_checkdata(freqPost,'cmbrepresentation', 'full');
else
    TfreqPost = freqPost;
end

cfg=[];
switch method
    case 'granger_bivariate'
        cfg.method = 'granger';
        cfg.granger.sfmethod = 'bivariate';
        cfg.granger.conditional = 'no';
        cfg.channelcmb = {'all','all'};
    case {'granger_multivariate','granger'}
        cfg.method = 'granger';
        cfg.granger.sfmethod = 'multivariate';
        cfg.granger.conditional = 'no';
    otherwise
        cfg.method = method;
end
fieldtripMethod = cfg.method;
conn_stat = ft_connectivityanalysis(cfg, TfreqPost);
end
