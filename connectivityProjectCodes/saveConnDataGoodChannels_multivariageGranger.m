function saveConnDataGoodChannels(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange,connDataFolder)
% Copy of saveConnData for multivariate GC using only good electrodes.
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

        [connPre,freqPre]=computeGC(dataPre,labels);
        [connPost,freqPost]=computeGC(dataPost,labels);
    end

    % Explicit provenance distinguishes these files from the original results.
    gcSettings=struct('sfmethod','multivariate','conditional','no', ...
        'badChannelsRemovedBeforeGC',true,'matrixOrder','source_target_freq', ...
        'stRange',stRange,'tapsmofrq',1,'foilim',[0 200]);
    outputFile=fullfile(subjectFolder, ...
        [protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);
    save(outputFile,'connPre','connPost','freqPre','freqPost','numGoodTrials', ...
        'labels','badElecs','goodElecs','gcSettings','inputFile');
    fprintf('Saved %s\n',outputFile);
end
end

function [conn,freqVals]=computeGC(data,labels)
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
cfg.granger.sfmethod='multivariate';
cfg.granger.conditional='no'; % Preserve original non-conditional GC setting.
G=ft_connectivityanalysis(cfg,freqData);

% Use returned labels, not assumed output order, to restore both matrix axes.
[~,idx]=ismember(G.label(:),labels);
freqVals=G.freq;
conn=nan(64,64,numel(freqVals));
conn(idx,idx,:)=G.grangerspctrm; % Keep raw GC; do not hide negatives using abs().
end
