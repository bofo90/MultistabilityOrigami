%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Matlab version used to write the code: R2014b
%
%FOR ANY QUESTION RELATED TO THE MATLAB FILES, PLEASE CONTACT JOHANNES
%OVERVELDE AT J.T.B.OVERVELDE@GMAIL.COM (WWW.OVERVELDE.COM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
format long
clearvars -global
addpath([pwd  '\Modules'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE PREDEFINED GEOMETRY, SIMULATION AND PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt=initOpt('inputType','individual', 'template','truncated tetrahedron', ...
            'analysis', 'plot', 'readHingeFile', 'off', 'Kappa',0.0001);    

opt.saveFile = strcat('/',date,'_Example');
opt.hingeSet = [1 2]; %This vector can be changed for the hinge numbers that you want to fold or plot

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findDeformation(unitCell,extrudedUnitCell,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReadAndPlot(unitCell, extrudedUnitCell, opt);

t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


