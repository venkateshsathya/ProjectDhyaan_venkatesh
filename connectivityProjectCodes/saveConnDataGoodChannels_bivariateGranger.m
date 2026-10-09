function saveConnDataGoodChannels_bivariate(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder)
% Bivariate counterpart of saveConnDataGoodChannels, using only good electrodes.
% Same calling arguments and saved connPre/connPost format as saveConnData.
% Use connMethod='granger' and a NEW connDataFolder to keep old results.
% Output: 64 x 64 x frequency, in the original data.label order.
% GC(source,target,f); excluded rows/columns remain NaN.
% Requires FieldTrip on the MATLAB path. Bad trials were removed upstream.

if nargin<7 || isempty(stRange), stRange=[0.25 1.25]; end
protocolNameList=cellstr(string(protocolNameList));
subjectFolder=fullfile(connDataFolder,subjectName);
if ~isfolder(subjectFolder), mkdir(subjectFolder); end

for p=1:numel(protocolNameList)
    protocolName=protocolNameList{p};
    inputFile=fullfile(ftDataFolder,subjectName, ...
        [protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
    S=load(inputFile);
    numGoodTrials=S.numGoodTrials;
    labels={}; badElecs=[]; goodElecs=[];
    connPre=[]; connPost=[]; freqPre=[]; freqPost=[];

    if numGoodTrials>0
        data=S.data;
        labels=data.label(:); % Original full 64-channel order, before selection.
        badElecs=data.badElecs;
        if islogical(badElecs), badElecs=find(badElecs); end
        badElecs=unique(badElecs(:));
        goodElecs=setdiff((1:64)',badElecs,'stable');

        % Remove custom bookkeeping fields that FieldTrip may treat as data.
        data=rmfield(data,intersect(fieldnames(data),{'badElecs','timeVals','numTrials'}));

        % Remove bad signals BEFORE Fourier analysis and spectral factorization.
        cfg=[];
        cfg.channel=labels(goodElecs);
        dataGood=ft_selectdata(cfg,data);

        % Preserve the original Pre/Post time windows.
        cfg=[];
        cfg.toilim=[-diff(stRange)+1/data.fsample 0];
        dataPre=ft_redefinetrial(cfg,dataGood);
        cfg.toilim=[stRange(1)+1/data.fsample stRange(2)];
        dataPost=ft_redefinetrial(cfg,dataGood);

        [connPre,freqPre]=computeGC(dataPre,labels,goodElecs);
        [connPost,freqPost]=computeGC(dataPost,labels,goodElecs);
    end

    % Explicit provenance distinguishes these files from the original results.
    gcSettings=struct('sfmethod','bivariate','conditional','no', ...
        'badChannelsRemovedBeforeGC',true,'matrixOrder','source_target_freq', ...
        'stRange',stRange,'tapsmofrq',1,'foilim',[0 200]);
    outputFile=fullfile(subjectFolder, ...
        [protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);
    save(outputFile,'connPre','connPost','freqPre','freqPost','numGoodTrials', ...
        'labels','badElecs','goodElecs','gcSettings','inputFile');
    fprintf('Saved %s\n',outputFile);
end
end

function [conn,freqVals]=computeGC(data,labels,goodElecs)
cfg=[];
cfg.method='mtmfft';
cfg.taper='dpss';
cfg.tapsmofrq=1;
cfg.output='fourier';
cfg.keeptrials='yes';
cfg.foilim=[0 200];
freqData=ft_freqanalysis(cfg,data);

cfg=[];
cfg.method='granger';
cfg.granger.sfmethod='bivariate';
cfg.granger.conditional='no';
cfg.channelcmb={'all','all'}; % FieldTrip handles all pairs in this ONE call.
G=ft_connectivityanalysis(cfg,freqData);

% Pair rows are source -> target; local FieldTrip labels include [pair] suffixes.
pairLabels=regexprep(G.labelcmb,'\[.*$','');
[~,sourceIdx]=ismember(pairLabels(:,1),labels);
[~,targetIdx]=ismember(pairLabels(:,2),labels);
freqVals=G.freq;
conn=nan(64,64,numel(freqVals));
for k=1:size(pairLabels,1)
    conn(sourceIdx(k),targetIdx(k),:)=reshape(G.grangerspctrm(k,:),1,1,[]);
end
% FieldTrip omits self-pairs. Match the original saver and seed-validity check.
for k=goodElecs(:).'
    conn(k,k,:)=0;
end
end
