function opt=initOpt(varargin)
%initOpt() initialize the simulations with default options. To change 
%specific options use initOpt('option1',value1,'option2',value2,...).
%
%The most important options are:
%1) 'template': any string indicating the specific template used to generate 
%   the architected material. 
%   -> When 'individual' is used, possible values are the names of the
%   platonic solids, archimedean solids and some prisms, including: 
%   'tetrahedron', 'cube', 'octahedron', truncated tetrahedron',
%   'cuboctahedron', 'truncated cube', 'truncated octahedron',
%   'rhombicuboctahedron', 'truncated cuboctahedron', 'triangular prism',
%   'hexagonal prism', 'octagonal prism' and 'decagonal prism'.
%2) 'analysis': 'selecthinges', 'info', 'result', 'savedata' or 'plot'. 
%   'selecthinges' runs the program to select all possible hinges for a
%   defined template and creates a file with this selection. 'info' should
%   be used when setting up a new template, while 'result' will run the 
%   simulation and fold the structure saving the results in .mat files. 
%   'savedata' will read the results, save the analysis on .csv files and 
%   can plot the defomrmation of the structure. 'plot' only does the latter.
%3) 'readHingeFile': 'on' or 'off'. This hinge file can be created by
%   selecting the analysis selecthinges. If it already exists, by turning 
%   this option 'on' all the possible hinge selections will be read from 
%   the file. Turn it 'off' and the user will have to provide a hinge 
%   selection. Works for simulating a deformation and for plotting it. 
%
%See the intiOpt() file for a few other (less important) options. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GEOMETRY OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type of geometry ('individual' or 'material')
opt.inputType = 'individual';
%Name of the template (name of polyhedron or number of material)
opt.template = 'cube';
%Apply periodic boundary conditions (on or off) only when material is
%selected
opt.periodic = 'on';
%Analysis type to do (info, result, savedata, plot)
opt.analysis = 'result';
%Extrusion length
%Value >0
opt.Lextrude=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stiffness of the different springs
opt.Khinge = 0.001;
opt.Kedge = 1;
opt.Kdiag = 1;
opt.Kface = 100;
opt.KtargetAngle = 100;

%Read file with hinge combinations to fold (on or off)
opt.readHingeFile = 'off';
%Specify hinges to fold and angle to fold (-pi is completely closed and 0 
%is open)
opt.hingeSet = [1 2];

%Specify the contraints on the structures
opt.constrEdge = 'off';
opt.constrFace = 'on';
opt.constAnglePerc = 0.99; 
opt.maxStretch = 0.3;

%Minimization algorithm
opt.folAlgor = 'sqp';
opt.relAlgor = opt.folAlgor;
%Fold structure in specified amount of steps
opt.steps = 1;
%Save intermidiate stepsduring minimization (on or off
opt.gethistory = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT AND MOVIE OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type of results to display         
%Number of frames to show the modes the modes
%Integer value > 0
opt.frames=30;
%Figure DPI
opt.figDPI = 200;

%When using periodic boundary condtions, number of unit cells to show
%Integer value >0
opt.plotPer=2;
%Plot options default view
opt.AZ=-54;
opt.EL=14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPDATE OPTIONS BASED ON INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2:nargin
    if strcmp(varargin{i},'Kappa')
        opt.Khinge = varargin{i+1};
    end
    opt.(varargin{i}) = varargin{i+1};
end
