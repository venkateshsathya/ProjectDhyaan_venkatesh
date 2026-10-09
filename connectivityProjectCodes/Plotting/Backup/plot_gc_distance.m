function [results, figs] = plot_gc_distance(dataRoot, subject, protocols, varargin)
%PLOT_GC_DISTANCE Saved ProjectDhyaan GC versus six cosine-distance bins.
% [R,F] = plot_gc_distance(root,'019CKa',{'EO1'},'Direction','both');
% root contains subject folders. No FieldTrip/EEGLAB toolbox is required.
% Verified saved schema: connPre/connPost [channel channel frequency],
% freqPre/freqPost, numGoodTrials. Companion FT data supplies labels/geometry.
% FieldTrip GC convention: C(source,target,f), verified in its implementation.
% This loader assumes the saver retained the companion data.label order;
% the inspected saver does so, but does not store labels in the GC MAT file.
%
% Options:
% FTRoot        sibling data/ftData by default
% Band          [20 32] Hz, inclusive; equally weighted sampled frequencies
% Epoch         auto | pre | post | mean. auto: post for G1/G2/M2,
%               mean of Pre and Post for EO1/EC1/M1/EO2/EC2 (paper convention)
% Direction     both | outgoing | incoming (never silently symmetrized)
% Distance      project | spherical. project reproduces the local helper's
%               wrapped azimuth/elevation norm and Cz adjustment. spherical
%               is the normalized 3-D dot product (true cos angular distance).
% AngleFile     label/angle CSV for project mode; bundled actiCap64_UOL
% MatrixOrder   source-target (FieldTrip default) | target-source
% Suffix        ep_v8_granger.mat
% FTSuffix      ep_v8.mat
% Visible       on | off
%
% Self-connections excluded. NaN/Inf are missing, never zero. A pair must
% have all selected frequency samples and requested epochs to contribute.
% Seeds with missing bins stay NaN. Black mean requires all four seed means.
% Returns matrices of per-seed/per-electrode values, memberships and counts.

p = inputParser;
p.addParameter('FTRoot',fullfile(fileparts(char(dataRoot)),'data','ftData'));
p.addParameter('Band',[20 32]);
p.addParameter('Epoch','auto');
p.addParameter('Direction','both');
p.addParameter('Distance','project');
p.addParameter('MatrixOrder','source-target');
p.addParameter('AngleFile',fullfile(fileparts(mfilename('fullpath')),'actiCap64_UOL_angles.csv'));
p.addParameter('Suffix','ep_v8_granger.mat');
p.addParameter('FTSuffix','ep_v8.mat');
p.addParameter('Visible','on');
p.parse(varargin{:}); o=p.Results;
o.Epoch=validatestring(o.Epoch,{'auto','pre','post','mean'});
o.Direction=validatestring(o.Direction,{'both','outgoing','incoming'});
o.Distance=validatestring(o.Distance,{'project','spherical'});
o.MatrixOrder=validatestring(o.MatrixOrder,{'source-target','target-source'});
validateattributes(o.Band,{'numeric'},{'vector','numel',2,'finite','nonnegative'});
assert(o.Band(2)>=o.Band(1),'Band must be increasing.');
protocols=cellstr(string(protocols)); subject=char(subject);
seeds={'Oz','O2','POz','O1'}; edges=linspace(-1,1,7);
centers=(edges(1:end-1)+edges(2:end))/2;
results=struct([]);
for q=1:numel(protocols)
    protocol=protocols{q};
    gcFile=fullfile(dataRoot,subject,[protocol '_' char(o.Suffix)]);
    ftFile=fullfile(o.FTRoot,subject,[protocol '_' char(o.FTSuffix)]);
    assert(isfile(gcFile),'Missing GC file: %s',gcFile);
    assert(isfile(ftFile),'Missing companion FT file with labels/geometry: %s',ftFile);
    S=load(gcFile); T=load(ftFile,'data');
    assert(isfield(T,'data') && isfield(T.data,'label') && isfield(T.data,'elec'), ...
        'Companion file must contain data.label and data.elec.');
    labels=cellstr(string(T.data.label(:))); n=numel(labels);
    assert(numel(unique(lower(string(labels))))==n,'Duplicate channel labels.');
    E=T.data.elec;
    assert(isfield(E,'label') && isfield(E,'chanpos'),'Missing elec.label/chanpos.');
    [found,ix]=ismember(lower(string(labels)),lower(string(E.label(:))));
    assert(all(found),'Some data.label entries have no electrode geometry.');
    xyz=double(E.chanpos(ix,:));
    assert(isequal(size(xyz),[n 3]) && all(isfinite(xyz(:))), 'Invalid electrode positions.');
    radii=sqrt(sum(xyz.^2,2)); assert(all(radii>0),'Zero-radius electrode coordinate.');
    unit=xyz./radii; % Origin is the cap sphere origin in these saved files.
    [found,seedIdx]=ismember(lower(string(seeds)),lower(string(labels)));
    assert(all(found),'Missing one or more seeds Oz, O2, POz, O1.');
    ep=o.Epoch;
    if strcmp(ep,'auto')
        if ismember(upper(protocol),{'G1','G2','M2'}), ep='post';
        elseif ismember(upper(protocol),{'EO1','EC1','M1','EO2','EC2'}), ep='mean';
        else, error('Unknown protocol: set Epoch explicitly.'); end
    end
    if strcmp(ep,'pre')
        [C,f]=readEpoch(S,'Pre',n,o.Band);
    elseif strcmp(ep,'post')
        [C,f]=readEpoch(S,'Post',n,o.Band);
    else
        [A,f]=readEpoch(S,'Pre',n,o.Band);
        [B,g]=readEpoch(S,'Post',n,o.Band);
        assert(isequal(f,g),'Pre/Post frequency grids differ.');
        C=(A+B)/2; % Strict: both epoch estimates required per pair/frequency.
    end
    if strcmp(o.MatrixOrder,'target-source'), C=permute(C,[2 1 3]); end
    C(~isfinite(C))=NaN;
    bandMean=mean(C,3); % No omitnan: identical frequency support per pair.
    bandMean(1:n+1:end)=NaN;
    if isfield(T.data,'badElecs')
        bad=T.data.badElecs;
        if islogical(bad), assert(numel(bad)==n,'Invalid badElecs mask.'); bad=find(bad); end
        assert(isnumeric(bad) && all(bad(:)>=1 & bad(:)<=n & mod(bad(:),1)==0), ...
            'Expected badElecs to contain channel indices.');
        bandMean(bad,:)=NaN; bandMean(:,bad)=NaN;
    end
    R=struct; R.protocol=protocol; R.subject=subject; R.epoch=ep;
    R.options=o; R.gcFile=gcFile; R.ftFile=ftFile;
    R.labels=labels; R.seeds=seeds; R.seedIndices=seedIdx;
    R.frequencies=f; R.edges=edges; R.centers=centers;
    R.cosDistance=nan(4,n); R.bin=nan(4,n);
    R.outgoing=nan(4,6); R.incoming=nan(4,6);
    R.countOutgoing=zeros(4,6); R.countIncoming=zeros(4,6);
    R.pairOutgoing=bandMean(seedIdx,:);
    R.pairIncoming=bandMean(:,seedIdx).';
    [az,el]=cart2sph(xyz(:,1),xyz(:,2),xyz(:,3));
    az=rad2deg(az); el=rad2deg(el);
    if strcmp(o.Distance,'project')
        assert(isfile(o.AngleFile),'Project mode requires AngleFile: %s',o.AngleFile);
        angles=readtable(o.AngleFile,'TextType','string');
        assert(numel(unique(lower(angles.Label)))==height(angles),'Duplicate angle labels.');
        [ok,mi]=ismember(lower(string(labels)),lower(angles.Label));
        assert(all(ok),'AngleFile lacks some channel labels; supply matching cap angles.');
        az=angles.AzimuthDegrees(mi); el=angles.ElevationDegrees(mi);
        assert(all(isfinite(az)) && all(isfinite(el)),'Invalid montage angles.');
    end
    for s=1:4
        i=seedIdx(s);
        if strcmp(o.Distance,'spherical')
            x=unit*unit(i,:).';
        else
            da=abs(mod(az-az(i)+180,360)-180);
            de=abs(mod(el-el(i)+180,360)-180);
            da(strcmpi(labels,'Cz'))=0; % Label-based replacement for hardcoded 24.
            x=cosd(hypot(da,de));
        end
        x=max(-1,min(1,x)); b=discretize(x,edges); b(i)=NaN;
        R.cosDistance(s,:)=x.'; R.bin(s,:)=b.';
        for k=1:6
            vo=R.pairOutgoing(s,b==k); vi=R.pairIncoming(s,b==k);
            R.countOutgoing(s,k)=sum(isfinite(vo));
            R.countIncoming(s,k)=sum(isfinite(vi));
            R.outgoing(s,k)=mean(vo,'omitnan');
            R.incoming(s,k)=mean(vi,'omitnan');
        end
    end
    R.meanOutgoing=mean(R.outgoing,1); R.meanIncoming=mean(R.incoming,1);
    if any(R.countOutgoing(:)==0) || any(R.countIncoming(:)==0)
        warning('plot_gc_distance:MissingBins', ...
            '%s: empty seed/bin estimates remain NaN; inspect countOutgoing/countIncoming.',protocol);
    end
    if isfield(S,'numGoodTrials'), R.numGoodTrials=S.numGoodTrials; else, R.numGoodTrials=NaN; end
    if q==1
        % Establish the field schema before indexed struct assignment.
        results=R;
    else
        results(q)=R; %#ok<AGROW>
    end
end
if strcmp(o.Direction,'both'), dirs={'outgoing','incoming'}; else, dirs={o.Direction}; end
figs=gobjects(1,numel(dirs));
colors=[0 .447 .741; .85 .325 .098; .929 .694 .125; .494 .184 .556];
for d=1:numel(dirs)
    dir=dirs{d};
    figs(d)=figure('Color','w','Visible',o.Visible,'Position',[70 60 1100 850]);
    if numel(protocols)==1, nr=1; nc=1; else, nc=2; nr=ceil(numel(protocols)/2); end
    tl=tiledlayout(figs(d),nr,nc,'TileSpacing','compact','Padding','compact');
    for q=1:numel(results)
        R=results(q); ax=nexttile(tl); hold(ax,'on'); Y=R.(dir);
        for s=1:4
            plot(ax,centers,Y(s,:),'-o','Color',colors(s,:), ...
                'LineWidth',1.3,'MarkerSize',4,'DisplayName',seeds{s});
        end
        plot(ax,centers,mean(Y,1),'-ok','LineWidth',2,'MarkerSize',5,'DisplayName','4-seed mean');
        set(ax,'XDir','reverse','XLim',[-1 1],'XTick',[-1 -.5 0 .5 1]);
        xlabel(ax,'cos\theta (near \rightarrow far)'); ylabel(ax,'Granger causality');
        title(ax,sprintf('%s | %s',R.protocol,R.epoch),'Interpreter','none');
        grid(ax,'on'); box(ax,'off');
        if ~any(isfinite(Y(:)))
            text(ax,.5,.5,'No valid seed-to-electrode estimates','Units','normalized', ...
                'HorizontalAlignment','center');
        end
        if q==1, legend(ax,'Location','best'); end
    end
    if strcmp(dir,'outgoing'), directionLabel='seed -> other electrodes';
    else, directionLabel='other electrodes -> seed'; end
    title(tl,sprintf('%s | %g-%g Hz | %s | %s distance',subject, ...
        o.Band(1),o.Band(2),directionLabel,o.Distance),'Interpreter','none');
end
end

function [C,f]=readEpoch(S,epoch,n,band)
cn=['conn' epoch]; fn=['freq' epoch];
assert(isfield(S,cn) && isfield(S,fn),'Missing fields %s / %s.',cn,fn);
C=S.(cn); f=double(S.(fn)(:).');
assert(isnumeric(C) && isreal(C) && ~isempty(C), ...
    '%s must be a nonempty real numeric channel-channel-frequency array.',cn);
assert(size(C,1)==n && size(C,2)==n && size(C,3)==numel(f) && ndims(C)<=3, ...
    '%s shape disagrees with companion labels or frequencies.',cn);
assert(all(isfinite(f)) && all(diff(f)>0),'Invalid frequency axis.');
assert(band(1)>=f(1) && band(2)<=f(end),'Band exceeds available frequencies.');
use=f>=band(1) & f<=band(2); assert(any(use),'No frequencies in band.');
C=double(C(:,:,use)); f=f(use);
end
