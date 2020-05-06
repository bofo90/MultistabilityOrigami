function metadataFile(opt, unitCell, extrudedUnitCell)

extraName = sprintf('/kh%2.3f_kta%2.3f_ke%2.3f_kf%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge, opt.Kface);
filedir = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat', opt.saveFile,extraName, '/');
filenametxt = strcat(filedir, 'metadata.txt');
filenamemat = strcat(filedir, 'metadata.mat');
opt.date = datestr(now(), 'yyyy-mm-dd HH:MM:SS');
[~, gitRev] = system('git rev-parse --verify --short HEAD');
opt.gitVersion = gitRev;
    
if strcmp(opt.analysis,'result')
    metadata.minimizationOpt = opt.options;
    metadata.options = rmfield(opt,{'options','angleConstrFinal'});
%     metadata.UnitCell.nodes = size(unitCell.Polyhedron.node,1);
%     metadata.UnitCell.edges = size(unitCell.Polyhedron.edge,1);
%     metadata.UnitCell.faces = size(unitCell.Polyhedron.face,1);
    metadata.extUnitCell.nodes = size(extrudedUnitCell.node,1);
    metadata.extUnitCell.edges = size(extrudedUnitCell.edge,1)-size(extrudedUnitCell.diagonals,2);
    metadata.extUnitCell.diag = size(extrudedUnitCell.diagonals,2);
    metadata.extUnitCell.faces = size(extrudedUnitCell.face,2);
    metadata.extUnitCell.hinges = size(extrudedUnitCell.nodeHingeEx,1);
    metadata.extUnitCell.intHinges = size(extrudedUnitCell.innerHinges,2);
    metadata.extUnitCell.maxStretch = extrudedUnitCell.maxStretch;
    metadata.extUnitCell.maxHingeFold = extrudedUnitCell.maxHingeFold;

    if ~exist(filedir, 'file')
        mkdir(filedir);
    end

    struct2ini(filenametxt, metadata)
    save(filenamemat, 'metadata');
end

function struct2ini(filename,Structure)
%==========================================================================
% Author:      Dirk Lohse ( dirklohse@web.de )
% Version:     0.1a
% Last change: 2008-11-13
%==========================================================================
% Modified by Agustin Iniguez 30-11-2017
%
% struct2ini converts a given structure into an ini-file.
% It's the opposite to Andriy Nych's ini2struct. Only 
% creating an ini-file is implemented. To modify an existing
% file load it with ini2struct.m from:
%       Andriy Nych ( nych.andriy@gmail.com )
% change the structure and write it with struct2ini.
%

% Open file, or create new file, for writing
% discard existing contents, if any.
fid = fopen(filename,'wt'); 

Sections = fieldnames(Structure);                     % returns the Sections

for i=1:length(Sections)
   Section = char(Sections(i));                       % convert to character
   
   fprintf(fid,'\n[%s]\n',Section);                       % output [Section]
   
   member_struct = Structure.(Section);               % returns members of Section
   if ~isempty(member_struct)                         % check if Section is empty
      member_names = fieldnames(member_struct);
      for j=1:length(member_names)
         member_name = char(member_names(j));
         member_value = Structure.(Section).(member_name);
         
         if isa(member_value, 'function_handle')
             member_value = func2str(member_value); 
         end
         if isnumeric(member_value)
             member_value = num2str(member_value);
         end
         if islogical(member_value)
             if member_value
                 member_value = 'True';
             else
                 member_value = 'False';
             end
         end
         
         fprintf(fid,'%s=%s\n',member_name,member_value); % output member name and value
         
      end % for-END (Members)
   end % if-END
end % for-END (Sections)

fclose(fid); % close file
