% Put this script and plot_gc_distance.m in the same folder; run this script.
addpath(fileparts(mfilename('fullpath')));
dataRoot='/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-08_17-43-56';
dataRoot = '/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-09_14-49-12';

subject='019CKa';
protocols={'EO1'}; % Selected protocol, with four seed curves + black mean.
% Uncomment for the 4 x 2 protocol grid (one figure for each GC direction):
% protocols={'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};

[results,figs]=plot_gc_distance(dataRoot,subject,protocols, ...
    'Band',[20 32],'Epoch','auto','Direction','both','Distance','project');

% Optional true spherical angular distance instead of the project's helper:
% [results,figs]=plot_gc_distance(dataRoot,subject,protocols,'Distance','spherical');

% Optional export (choose an output folder):
% outDir=fullfile(fileparts(mfilename('fullpath')),'plots');
% if ~isfolder(outDir), mkdir(outDir); end
% exportgraphics(figs(1),fullfile(outDir,[subject '_outgoing.png']),'Resolution',200);
% exportgraphics(figs(2),fullfile(outDir,[subject '_incoming.png']),'Resolution',200);
% save(fullfile(outDir,[subject '_binned_gc.mat']),'results');
