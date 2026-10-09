function test_gc_distance
% Synthetic directed fixture: run test_gc_distance in MATLAB.
% Checks orientation, label-based geometry, exact endpoints and self removal.
r=tempname; mkdir(r); cleanup=onCleanup(@()rmdir(r,'s')); %#ok<NASGU>
mkdir(fullfile(r,'S')); mkdir(fullfile(r,'ft','S'));
labels={'Oz';'O2';'POz';'O1';'Near';'Far'};
data.label=labels;
xyz=[1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;-1 0 0];
perm=[6 4 2 5 1 3]; data.elec.label=labels(perm); data.elec.chanpos=xyz(perm,:);
data.badElecs=[];
save(fullfile(r,'ft','S','EO1_ep_v8.mat'),'data');
connPre=repmat((1:6)'*10+(1:6),1,1,2); connPost=connPre+2;
freqPre=[20 32]; freqPost=freqPre; numGoodTrials=10;
save(fullfile(r,'S','EO1_ep_v8_granger.mat'),'connPre','connPost','freqPre','freqPost','numGoodTrials');
[R,F]=plot_gc_distance(r,'S','EO1','FTRoot',fullfile(r,'ft'),'Distance','spherical','Visible','off');
close(F);
assert(R.countOutgoing(1,6)==4 && R.countOutgoing(1,1)==1);
assert(isnan(R.bin(1,1)) && R.bin(1,5)==6 && R.bin(1,6)==1);
assert(R.outgoing(1,1)==17 && R.incoming(1,1)==62);
assert(R.outgoing(1,6)==14.5 && R.incoming(1,6)==37);
assert(all(isnan(R.outgoing(:,2:5)),'all'));
% Regression: first-result initialization and subsequent protocol assignment.
copyfile(fullfile(r,'ft','S','EO1_ep_v8.mat'),fullfile(r,'ft','S','EC1_ep_v8.mat'));
save(fullfile(r,'S','EC1_ep_v8_granger.mat'),'connPre','connPost','freqPre','freqPost');
[Q,F]=plot_gc_distance(r,'S',{'EO1','EC1'},'FTRoot',fullfile(r,'ft'), ...
    'Distance','spherical','Visible','off'); close(F);
assert(numel(Q)==2 && strcmp(Q(1).protocol,'EO1') && strcmp(Q(2).protocol,'EC1'));
assert(isequaln(Q(1).outgoing,Q(2).outgoing));
assert(Q(1).numGoodTrials==10 && isnan(Q(2).numGoodTrials));
% Project-mode CSV is also matched by labels, not row number.
Label=labels(perm); AzimuthDegrees=[0;0;0;0;0;180];
AzimuthDegrees=AzimuthDegrees(perm); ElevationDegrees=zeros(6,1);
angleFile=fullfile(r,'angles.csv');
writetable(table(Label,AzimuthDegrees,ElevationDegrees),angleFile);
[P,F]=plot_gc_distance(r,'S','EO1','FTRoot',fullfile(r,'ft'), ...
    'Distance','project','AngleFile',angleFile,'Visible','off'); close(F);
assert(isequaln(P.outgoing,R.outgoing) && isequaln(P.bin,R.bin));
% The same physical GC stored with target as first dimension must agree.
connPre=permute(connPre,[2 1 3]); connPost=permute(connPost,[2 1 3]);
save(fullfile(r,'S','EO1_ep_v8_granger.mat'),'connPre','connPost','freqPre','freqPost');
[B,F]=plot_gc_distance(r,'S','EO1','FTRoot',fullfile(r,'ft'), ...
    'MatrixOrder','target-source','Distance','spherical','Visible','off'); close(F);
assert(isequaln(R.outgoing,B.outgoing) && isequaln(R.incoming,B.incoming));
% A missing sample must invalidate that pair, not bias its band average.
connPre(6,1,1)=NaN;
save(fullfile(r,'S','EO1_ep_v8_granger.mat'),'connPre','connPost','freqPre','freqPost');
[B,F]=plot_gc_distance(r,'S','EO1','FTRoot',fullfile(r,'ft'), ...
    'MatrixOrder','target-source','Distance','spherical','Visible','off'); close(F);
assert(isnan(B.outgoing(1,1)) && B.countOutgoing(1,1)==0);
fprintf('All directed GC, geometry, bin-edge and missing-data tests passed.\n');
end
