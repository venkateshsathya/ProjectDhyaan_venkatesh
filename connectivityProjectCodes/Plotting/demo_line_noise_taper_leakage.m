% DEMO_LINE_NOISE_TAPER_LEAKAGE
% Reproducible PSD demonstration using the settings in saveConnData.m:
% mtmfft, DPSS, tapsmofrq=1 Hz, 1-second epochs, native 1-Hz FFT bins.
% Compare a rectangular window (no taper) with the pipeline DPSS taper.
% The line-only signal isolates leakage without clean/line cross terms.
% No connectivity or Granger causality is estimated by this demonstration.
%
% Run this script from MATLAB. Requires dpss (Signal Processing Toolbox) and
% FieldTrip on the path, as in the connectivity pipeline.

%% Inputs
fs = 1000;                  % Hz; assumed sampling rate, change to your data.fsample.
stRange = [0.25 1.25];       % Same default as saveConnData.m; duration is 1 second.
tapsmofrq = 1;              % Hz, DPSS half-bandwidth (pipeline value).
numTrials = 200;
randomSeed = 27;
cleanRmsUv = 5;             % RMS of the simulated clean signal, in microvolts.
lineToCleanPowerDb = 20;    % +20 dB = 100x power = 10x RMS amplitude.
lineFrequenciesHz = [50 50.3]; % On-bin line and small frequency offset.
strengthSweepDb = -20:5:40;
lineRemovalHz = [50 100];   % Same fitted-sinusoid removal as the DFT pipeline.
analysisMaxHz = 200;        % Same displayed/analysed range as saveConnData.m.
plotMaxHz = 120;
verifyAgainstFieldTrip = true; % Check the explicit PSD calculation against mtmfft.
saveOutputs = true;
outputFolder = fullfile(fileparts(mfilename('fullpath')),'line_noise_taper_demo_results');

%% Exact taper construction and PSD convention
% FieldTrip calls dpss(N, N*tapsmofrq/fs), then drops its final sequence.
% T=1 s and W=1 Hz -> NW=1 and K=1 retained taper (not several tapers).
% For each unit-energy taper w, one-sided PSD at interior frequencies is
%   P(f_k) = 2/(fs*sum(w.^2)) * |FFT((x-mean(x)).*w)|^2.
% DC and Nyquist are not doubled. Average over tapers and then trials.
% Native-grid band power = sum(P(selected bins))*df, in microvolt^2.
% A PSD integral measures mean-square power; multiply by T for epoch energy.
% Multiplication by w convolves the spectrum with W: nearby broadening is
% expected, but distant leakage need not increase. Compare both scenarios.
% Sources:
% https://www.fieldtriptoolbox.org/tutorial/sensor/timefrequencyanalysis/
% https://www.mathworks.com/help/signal/ref/dpss.html

assert(exist('dpss','file') ~= 0,'dpss is required (Signal Processing Toolbox).');
assert(exist('ft_preproc_dftfilter','file') ~= 0,'Add FieldTrip to the MATLAB path first.');
validateattributes(fs,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(numTrials,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cleanRmsUv,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(lineToCleanPowerDb,{'numeric'},{'scalar','real','finite'});
validateattributes(lineFrequenciesHz,{'numeric'},{'vector','real','finite','positive','<',fs/2});
validateattributes(tapsmofrq,{'numeric'},{'scalar','real','finite','positive'});
T = diff(stRange);
assert(isscalar(T) && T > 0,'stRange must contain an increasing pair of times.');
N = round(T*fs);
assert(abs(N/fs-T) < 1e-10 && mod(N,2)==0,'Use an even integer number of samples per epoch.');
assert(analysisMaxHz <= fs/2 && all(lineRemovalHz < fs/2),'Analysis/removal frequencies must fit below Nyquist.');
t = (0:N-1)/fs;
df = fs/N;
NW = N*tapsmofrq/fs;
[dpssTapers,concentration] = dpss(N,NW);
dpssTapers = dpssTapers(:,1:end-1);
concentration = concentration(1:end-1);
K = size(dpssTapers,2);
assert(K >= 1,'Epoch too short for this smoothing bandwidth.');
boxTaper = ones(N,1)/sqrt(N);
windows = {boxTaper,dpssTapers};
methodNames = {'Rectangular (no taper)','DPSS (pipeline)'};
freq = (0:N/2)*df;
analysisMask = freq > 0 & freq <= analysisMaxHz;
retainedMask = analysisMask & ~(freq >= 46 & freq <= 54);
[~,centerBin] = min(abs(freq-50));
outsideCenter = true(size(freq)); outsideCenter(centerBin) = false;
outsideSmoothing = abs(freq-50) > tapsmofrq;
outsideOmitted = freq < 46 | freq > 54;
fprintf('Epoch %.3f s; fs %.1f Hz; df %.3f Hz; NW %.3f; retained DPSS tapers %d.\n',T,fs,df,NW,K);

%% Simulate one clean dataset, reused for every noise condition
previousRng = rng;
rng(randomSeed,'twister');
burnIn = round(fs);
background = filter(1,[1 -0.95],randn(numTrials,N+burnIn),[],2);
background = background(:,burnIn+1:end);
background = background./std(background,0,2);
clean = background + 0.4*randn(numTrials,N) ...
    + 1.2*sin(2*pi*10*t+2*pi*rand(numTrials,1)) ...
    + 0.5*sin(2*pi*40*t+2*pi*rand(numTrials,1)) ...
    + 0.3*sin(2*pi*55*t+2*pi*rand(numTrials,1));
clean = clean-mean(clean,2);
clean = clean * cleanRmsUv/sqrt(mean(clean.^2,'all'));
% Evenly distributed phases avoid favoring one phase of the off-bin tone.
linePhases = 2*pi*randperm(numTrials)'/numTrials;
rng(previousRng);
cleanPower = mean(clean.^2,'all');
linePowerTarget = cleanPower*10^(lineToCleanPowerDb/10);
lineRmsUv = sqrt(linePowerTarget);
cleanFiltered = ft_preproc_dftfilter(clean,fs,lineRemovalHz,'dftreplace','zero');
cleanPSD = cellfun(@(w) estimatePSD(clean,w,fs,N),windows,'UniformOutput',false);
cleanFilteredPSD = estimatePSD(cleanFiltered,dpssTapers,fs,N);
% Parseval validation for the rectangular estimator, after demeaning.
assert(abs(sum(cleanPSD{1})*df-cleanPower) < 1e-10*max(1,cleanPower));
settings = struct('fs',fs,'stRange',stRange,'duration',T,'N',N,'df',df, ...
    'method','mtmfft','taper','dpss','tapsmofrq',tapsmofrq,'NW',NW,'numTapers',K, ...
    'numTrials',numTrials,'randomSeed',randomSeed,'cleanRmsUv',cleanRmsUv, ...
    'lineToCleanPowerDb',lineToCleanPowerDb,'lineFrequenciesHz',lineFrequenciesHz, ...
    'lineRemovalHz',lineRemovalHz,'analysisMaxHz',analysisMaxHz, ...
    'strengthSweepDb',strengthSweepDb,'lineRmsUv',lineRmsUv);
results = struct;
results.settings = settings;
results.freq = freq;
results.dpssTapers = dpssTapers;
results.concentration = concentration;
results.cleanPSD = cleanPSD;
results.cleanFilteredPSD = cleanFilteredPSD;
results.fieldTripRelativeError = NaN;
summaryRows = cell(0,9);
bandRows = cell(0,7);
bands = [8 12;30 45;46 50;51 55;55 80;80 120];
figs = gobjects(1,numel(lineFrequenciesHz)+1);
if saveOutputs && ~isfolder(outputFolder), mkdir(outputFolder); end
fprintf('Clean RMS = %.3f uV; line RMS = %.3f uV; power ratio = %.1f dB.\n', ...
    sqrt(cleanPower),lineRmsUv,lineToCleanPowerDb);

%% Matched comparisons, with and without contamination
for si = 1:numel(lineFrequenciesHz)
    fLine = lineFrequenciesHz(si);
    unitLine = sin(2*pi*fLine*t+linePhases);
    unitLine = unitLine-mean(unitLine,2);
    % Calibrate each trial to exactly unit mean-square power, including off-bin tones.
    unitLine = unitLine./sqrt(mean(unitLine.^2,2));
    lineOnly = lineRmsUv*unitLine;
    noisy = clean+lineOnly;
    filtered = ft_preproc_dftfilter(noisy,fs,lineRemovalHz,'dftreplace','zero');
    filteredLine = ft_preproc_dftfilter(lineOnly,fs,lineRemovalHz,'dftreplace','zero');
    if any(abs(fLine-lineRemovalHz) < 1e-12) && abs(fLine*T-round(fLine*T)) < 1e-12
        assert(norm(filtered-cleanFiltered,'fro')/norm(lineOnly,'fro') < 1e-10, ...
            'Known stationary line was not removed to numerical precision.');
    end
    assert(abs(mean(lineOnly.^2,'all')/cleanPower-10^(lineToCleanPowerDb/10)) < 1e-8);
    scenario = struct('lineHz',fLine,'lineRmsUv',sqrt(mean(lineOnly.^2,'all')), ...
        'lineToCleanPowerDb',10*log10(mean(lineOnly.^2,'all')/cleanPower));
    for mi = 1:2
        scenario.noisyPSD{mi} = estimatePSD(noisy,windows{mi},fs,N);
        scenario.linePSD{mi} = estimatePSD(lineOnly,windows{mi},fs,N);
    end
    scenario.filteredPSD = estimatePSD(filtered,dpssTapers,fs,N);
    scenario.filteredLinePSD = estimatePSD(filteredLine,dpssTapers,fs,N);
    % Line PSD is an exact isolated contribution. noisyPSD-cleanPSD also has
    % finite-sample cross terms and should not be interpreted as pure leakage.
    allLinePSD = [scenario.linePSD,{scenario.filteredLinePSD}];
    allNames = [methodNames,{'DFT removal then DPSS'}];
    references = [cleanPSD,{cleanFilteredPSD}];
    sweepScale = 10.^((strengthSweepDb-lineToCleanPowerDb)/10);
    scenario.sweepRetainedDb = zeros(3,numel(strengthSweepDb));
    for mi = 1:3
        P = allLinePSD{mi}; total = sum(P)*df;
        originalTotal = sum(scenario.linePSD{min(mi,2)})*df;
        retained = sum(P(retainedMask))*df;
        cleanRetained = sum(references{mi}(retainedMask))*df;
        % Fractions for filtered data use the unfiltered DPSS line power as denominator.
        summaryRows(end+1,:) = {fLine,allNames{mi},scenario.lineToCleanPowerDb,total, ...
            100*sum(P(outsideCenter))*df/originalTotal, ...
            100*sum(P(outsideSmoothing))*df/originalTotal, ...
            100*sum(P(outsideOmitted))*df/originalTotal,retained, ...
            100*retained/cleanRetained}; %#ok<SAGROW>
        scenario.sweepRetainedDb(mi,:) = 10*log10(max(retained*sweepScale/cleanRetained,realmin));
    end
    scenario.bandExcessPercent = zeros(size(bands,1),3);
    scenario.bandPowerRatioDb = zeros(size(bands,1),3);
    estimates = [scenario.noisyPSD,{scenario.filteredPSD}];
    % Performance relative to the SAME estimator applied to unfiltered clean data.
    % The filtered result includes any loss of clean 50/100 Hz signal power.
    cleanReferences = [cleanPSD,{cleanPSD{2}}];
    for bi = 1:size(bands,1)
        mask = freq >= bands(bi,1) & freq < bands(bi,2);
        for mi = 1:3
            cleanBand = sum(cleanReferences{mi}(mask))*df;
            observedBand = sum(estimates{mi}(mask))*df;
            bias = 100*(observedBand-cleanBand)/cleanBand;
            ratioDb = 10*log10(observedBand/cleanBand);
            scenario.bandExcessPercent(bi,mi) = bias;
            scenario.bandPowerRatioDb(bi,mi) = ratioDb;
            bandRows(end+1,:) = {fLine,allNames{mi},bands(bi,1),bands(bi,2),cleanBand,bias,ratioDb}; %#ok<SAGROW>
        end
    end
    % Validate native-grid DPSS PSD against the actual pipeline library.
    if verifyAgainstFieldTrip && si == 1
        assert(exist('ft_freqanalysis','file') ~= 0,'ft_freqanalysis is required for verification.');
        raw = struct('label',{{'simulated'}},'fsample',fs);
        nCheck = min(8,numTrials);
        raw.trial = mat2cell(noisy(1:nCheck,:),ones(1,nCheck),N).';
        raw.time = repmat({t},1,nCheck);
        raw.sampleinfo = [(0:nCheck-1)'*N+1,(1:nCheck)'*N];
        cfg = [];
        cfg.method = 'mtmfft'; cfg.taper = 'dpss'; cfg.tapsmofrq = tapsmofrq;
        cfg.output = 'pow'; cfg.foilim = [0 analysisMaxHz]; cfg.feedback = 'no';
        cfg.pad = T; cfg.polyremoval = 0; % Explicit versions of the pipeline defaults.
        ftResult = ft_freqanalysis(cfg,raw);
        explicitPSD = estimatePSD(noisy(1:nCheck,:),dpssTapers,fs,N);
        idx = round(ftResult.freq/df)+1;
        % FT unit-energy tapers: powspctrm is bin power; divide by df for PSD.
        ftPSD = reshape(ftResult.powspctrm,1,[])/df;
        relativeError = norm(ftPSD-explicitPSD(idx))/norm(ftPSD);
        assert(relativeError < 1e-9,'Explicit PSD does not match the active FieldTrip implementation.');
        results.fieldTripRelativeError = relativeError;
        results.fieldTripPath = which('ft_freqanalysis');
        fprintf('FieldTrip PSD cross-check relative error: %.3g\n',relativeError);
    end
    results.scenarios(si) = scenario;

    %% One figure for each line-frequency scenario
    figs(si) = figure('Color','w','Name',sprintf('Line noise %.1f Hz',fLine), ...
        'Position',[40 40 1500 900]);
    tl = tiledlayout(figs(si),2,3,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl);
    plot(ax,t,clean(1,:),'Color',[0 0.45 0.7]); hold(ax,'on');
    plot(ax,t,noisy(1,:),'Color',[0.85 0.3 0.1]); xlim(ax,[0 0.2]); grid(ax,'on');
    xlabel(ax,'Time (s)'); ylabel(ax,'Amplitude (uV)');
    title(ax,'Identical clean signal + calibrated line'); legend(ax,{'Clean','Contaminated'});

    ax = nexttile(tl); hold(ax,'on');
    semilogy(ax,freq,max(cleanPSD{1},realmin),'--','Color',[0 0.45 0.7]);
    semilogy(ax,freq,max(cleanPSD{2},realmin),'-','Color',[0 0.45 0.7]);
    semilogy(ax,freq,max(scenario.noisyPSD{1},realmin),'--','Color',[0.85 0.3 0.1]);
    semilogy(ax,freq,max(scenario.noisyPSD{2},realmin),'-','Color',[0.85 0.3 0.1]);
    set(ax,'YScale','log'); xlim(ax,[1 plotMaxHz]); grid(ax,'on');
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'PSD (uV^2/Hz)'); title(ax,'Clean/noisy x rectangular/DPSS');
    legend(ax,{'Clean rectangular','Clean DPSS','Noisy rectangular','Noisy DPSS'},'Location','best');

    ax = nexttile(tl); hold(ax,'on');
    % Floor only for display; leakage metrics always use unfloored PSD values.
    displayFloor = linePowerTarget*1e-12/df;
    semilogy(ax,freq,max(scenario.linePSD{1},displayFloor),'o--');
    semilogy(ax,freq,max(scenario.linePSD{2},displayFloor),'o-');
    set(ax,'YScale','log'); xlim(ax,[40 60]); grid(ax,'on');
    xline(ax,49,':'); xline(ax,51,':');
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Line-only PSD (uV^2/Hz)');
    title(ax,sprintf('Pure line: native %.2g-Hz bins',df));
    legend(ax,methodNames,'Location','best');

    ax = nexttile(tl); hold(ax,'on');
    semilogy(ax,freq,max(cleanPSD{2},realmin),'k-','LineWidth',1.4);
    semilogy(ax,freq,max(scenario.noisyPSD{2},realmin),'-','Color',[0.85 0.3 0.1]);
    semilogy(ax,freq,max(scenario.filteredPSD,realmin),'-','Color',[0 0.6 0.3]);
    semilogy(ax,freq,max(cleanFilteredPSD,realmin),'k:');
    set(ax,'YScale','log'); xlim(ax,[40 60]); grid(ax,'on');
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'PSD (uV^2/Hz)');
    title(ax,'Fitted 50/100 Hz removal BEFORE DPSS');
    legend(ax,{'Clean DPSS','Noisy DPSS','Filtered noisy DPSS','Filtered clean DPSS'},'Location','best');

    ax = nexttile(tl);
    leakagePercent = zeros(3,3);
    masks = {outsideCenter,outsideSmoothing,outsideOmitted};
    for mi = 1:3
        denominator = sum(scenario.linePSD{min(mi,2)})*df;
        for bi = 1:3
            leakagePercent(bi,mi) = 100*sum(allLinePSD{mi}(masks{bi}))*df/denominator;
        end
    end
    bar(ax,max(leakagePercent,1e-10)); set(ax,'YScale','log'); grid(ax,'on');
    xticks(ax,1:3); xticklabels(ax,{'Outside 50-Hz bin','Outside 50 +/- W','Outside 46-54 Hz'});
    ylabel(ax,'% of unfiltered line spectral power');
    title(ax,'Nested masks; full 0-Nyquist spectrum');
    legend(ax,allNames,'Location','best');

    ax = nexttile(tl);
    bar(ax,scenario.bandPowerRatioDb); grid(ax,'on');
    xticks(ax,1:size(bands,1)); xticklabels(ax,compose('%g-%g',bands(:,1),bands(:,2)));
    xlabel(ax,'Frequency band [low, high), Hz'); ylabel(ax,'Band power / clean reference (dB)');
    title(ax,'0 dB = clean; +10 dB = 10x power'); legend(ax,allNames,'Location','best');
    title(tl,sprintf('%.1f-Hz line | %.0f dB power ratio | T=%.1f s, DPSS W=%.1f Hz, K=%d', ...
        fLine,lineToCleanPowerDb,T,tapsmofrq,K));
    if saveOutputs
        drawnow; % Finish layout before export, including when figures are hidden.
        exportgraphics(tl,fullfile(outputFolder,sprintf('scenario_%gHz.png',fLine)),'Resolution',150);
    end
end

%% Taper shapes, spectral windows, and strength calibration
figs(end) = figure('Color','w','Name','Taper kernels and noise calibration','Position',[60 60 1500 480]);
tl = tiledlayout(figs(end),1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); plot(ax,t,boxTaper,'--'); hold(ax,'on'); plot(ax,t,dpssTapers);
xlabel(ax,'Time (s)'); ylabel(ax,'Unit-energy taper'); grid(ax,'on');
title(ax,sprintf('Pipeline retains %d DPSS taper(s)',K));
legend(ax,[{'Rectangular'},cellstr(compose('DPSS %d',1:K))],'Location','best');
ax = nexttile(tl); hold(ax,'on');
fineN = 32*N; offset = (-fineN/2:fineN/2-1)*fs/fineN;
for mi = 1:2
    kernel = mean(abs(fftshift(fft(windows{mi},fineN,1),1)).^2,2)/fs;
    plot(ax,offset,10*log10(max(kernel/max(kernel),1e-12)));
end
xlim(ax,[-10 10]); ylim(ax,[-100 2]); grid(ax,'on');
xline(ax,-tapsmofrq,':'); xline(ax,tapsmofrq,':');
xlabel(ax,'Offset from tone (Hz)'); ylabel(ax,'Spectral window relative to peak (dB)');
title(ax,'Dense kernel display; no added resolution'); legend(ax,methodNames);
ax = nexttile(tl); hold(ax,'on');
styles = {'--','-',':'}; scenarioColors = lines(numel(lineFrequenciesHz));
for si = 1:numel(lineFrequenciesHz)
    for mi = 1:3
        plot(ax,strengthSweepDb,results.scenarios(si).sweepRetainedDb(mi,:), ...
            'Color',scenarioColors(si,:),'LineStyle',styles{mi},'LineWidth',1.4, ...
            'DisplayName',sprintf('%g Hz: %s',lineFrequenciesHz(si),allNames{mi}));
    end
end
yline(ax,0,'k:','Line contribution equals clean power','HandleVisibility','off');
xlabel(ax,'Injected line / clean total power (dB)'); ylabel(ax,'Retained line / clean band power (dB)');
title(ax,'0-200 Hz excluding 46-54 Hz'); grid(ax,'on'); legend(ax,'show','Location','best','FontSize',8);
results.leakageTable = cell2table(summaryRows,'VariableNames', ...
    {'LineHz','Method','InjectedPowerDb','TotalLineSpectralPowerUv2','Outside50BinPercent', ...
    'OutsideSmoothingPercent','Outside46to54Percent','RetainedLinePowerUv2','RetainedLineVsCleanPercent'});
results.bandErrorTable = cell2table(bandRows,'VariableNames', ...
    {'LineHz','Method','BandLowHz','BandHighHz','CleanBandPowerUv2','BandPowerErrorPercent','BandPowerRatioDb'});
disp(results.leakageTable);
fprintf(['Fractions quantify isolated line power, not total contaminated-signal power.\n' ...
    'A 50-Hz fit need not remove a 50.3-Hz tone; inspect its residual.\n' ...
    'Tapering redistributes spectral power; it is not line-noise removal.\n']);
if saveOutputs
    drawnow;
    exportgraphics(tl,fullfile(outputFolder,'taper_kernels_and_strength_sweep.png'),'Resolution',150);
    writetable(results.leakageTable,fullfile(outputFolder,'leakage_summary.csv'));
    writetable(results.bandErrorTable,fullfile(outputFolder,'band_power_error.csv'));
    save(fullfile(outputFolder,'simulation_results.mat'),'results');
    fprintf('Saved results to %s\n',outputFolder);
end

function P = estimatePSD(x,tapers,fs,nfft)
% x = trials x samples; tapers = samples x K. Standard one-sided PSD.
x = x-mean(x,2);
P = zeros(1,nfft/2+1);
for k = 1:size(tapers,2)
    w = tapers(:,k).';
    X = fft(x.*w,nfft,2);
    pk = mean(abs(X(:,1:nfft/2+1)).^2,1)/(fs*sum(w.^2));
    pk(2:end-1) = 2*pk(2:end-1);
    P = P+pk/size(tapers,2);
end
end
