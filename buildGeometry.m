function [unitCell,extrudedUnitCell,opt]=buildGeometry(opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE EXTRUDED UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Get initial geometry
[unitCell,opt]=selectUnitCell(opt);
%Expand unit cell in case neccesary
[unitCell]=expandUnitCell(unitCell,opt);
%Create decoupled extruded faces
[extrudedUnitCell,unitCell]=extrudeUnitCell(unitCell,opt);

function [unitCell,opt]=selectUnitCell(opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch opt.inputType
    case {'user'}     
        %Initialize empty arrays
        unitCell.expCon=[];
        unitCell.perCon=[];
        %Load user file
        eval(opt.template)
        %Check for correct input, otherwise give error or warning message
        if isfield(unitCell,'Polyhedron')==0
            error('\n----------\nNo polyhedron specified for the template, make sure the userfile contains at least one call for each i-th polyhedron:\n\nPolyhedron(i)=polyhedra(''namePolyhedron'')\n----------\n',[]);
        end
        if and(length(unitCell.Polyhedron)>1,length(unitCell.expCon)+1<length(unitCell.Polyhedron))
            error(']n----------\nWhen the unit cell containts more than one polyhedron, for each i-th polyhedron (i=2:P, P being the total number of polyhedra in the unit cell) specify the location by giving congruent vertices between the different polyhedra:\n\nexpCon(i-1).dir=[v1a v1b p1;\n                 v2a v2b p2;\n                 v3a v3b p3].\n\nHere, each line represents one congruent vertex pair, in which v1a refers to the 1st vertex label of the i-th polyhedron, and v1b refers to the 1st vertex label of the p1-th polyhedron etc. Depending on the initial orientation of the polyhedra, the number of rows can be less than three.\n----------\n',[])
        end
        if size(unitCell.perCon,1)==0
            warning('\n----------\nNo periodic boundary conditions specified, running simulation without periodic constraints.\n\n To specify periodic boundary contions use:\n\nperCon=[f1a f1b p1a p1b;\n        f2a f2b p2a p2b;\n        f3a f3b p2a p2b].\n\nHere, f1a and f2a are the face pairs periodically located, which respectively belong to polyhedra with index p1a and p1b etc.\n----------\n',[])
        elseif size(unitCell.perCon,2)==2
            if length(unitCell.Polyhedron)>1
                warning('\n----------\nFor the face pair, no specific polyhedra are specified. Using the polyhedron with index 1 as default.\n----------\n',[])
            end
        end    
    case {'individual'}
        %Initialize empty array not used
        unitCell.perCon=[];
        unitCell.expCon=[];
        %Load polyhedron
        unitCell.Polyhedron(1)=polyhedra(opt.template);
    otherwise
        %LOAD ONE OF THE PREDEFINED UNIFORM SPACE-FILLING TESSELATIONS, OR 
        switch opt.template 
            case {'#1'}
                unitCell.Polyhedron(1)=polyhedra('tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.Polyhedron(3)=polyhedra('tetrahedron');
                unitCell.expCon(1).dir=[2 4 1; 1 2 1; 4 1 1];
                unitCell.expCon(2).dir=[3 6 2; 4 3 2; 2 5 2];
                unitCell.perCon=[1 1 1 2; 2 8 1 2; 4 5 1 2];
            case {'#2'}
                unitCell.Polyhedron(1)=polyhedra('tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.Polyhedron(3)=polyhedra('tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('octahedron');
                unitCell.Polyhedron(5)=polyhedra('tetrahedron');
                unitCell.Polyhedron(6)=polyhedra('tetrahedron');
                unitCell.expCon(1).dir=[2 4 1; 1 2 1; 4 1 1];
                unitCell.expCon(2).dir=[3 6 2; 4 3 2; 2 5 2];
                unitCell.expCon(3).dir=[4 4 2; 1 5 2; 2 6 2];
                unitCell.expCon(4).dir=[2 2 3; 4 1 3; 3 3 3];
                unitCell.expCon(5).dir=[2 5 4; 4 6 4; 3 4 4];
                unitCell.perCon=[1 3 1 6; 2 8 1 2; 4 5 1 2];          
            case {'#3'}
                unitCell.Polyhedron(1)=polyhedra('tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.Polyhedron(3)=polyhedra('tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.Polyhedron(5)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[2 4 1; 1 2 1; 4 1 1];
                unitCell.expCon(2).dir=[3 6 2; 4 3 2; 2 5 2];
                unitCell.expCon(3).dir=[5 4 2; 4 5 2; 6 6 2];
                unitCell.expCon(4).dir=[1 2 3; 3 3 3; 2 1 3];
                unitCell.perCon=[1 4 1 4; 2 8 1 2; 4 5 1 2];
            case {'#4'}
                unitCell.Polyhedron(1)=polyhedra('tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.Polyhedron(3)=polyhedra('tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.Polyhedron(5)=polyhedra('triangular prism');
                unitCell.Polyhedron(6)=polyhedra('octahedron');
                unitCell.Polyhedron(7)=polyhedra('tetrahedron');
                unitCell.Polyhedron(8)=polyhedra('tetrahedron');
                unitCell.Polyhedron(9)=polyhedra('triangular prism');
                unitCell.Polyhedron(10)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[2 4 1; 1 2 1; 4 1 1];
                unitCell.expCon(2).dir=[3 6 2; 4 3 2; 2 5 2];
                unitCell.expCon(3).dir=[5 4 2; 4 5 2; 6 6 2];
                unitCell.expCon(4).dir=[1 2 3; 3 3 3; 2 1 3];                    
                unitCell.expCon(5).dir=[5 1 4; 3 3 4; 1 2 4];                   
                unitCell.expCon(6).dir=[1 2 6; 4 4 6; 3 1 6];                   
                unitCell.expCon(7).dir=[2 4 5; 4 6 5; 1 5 5];
                unitCell.expCon(8).dir=[5 2 7; 4 4 7; 6 1 7];                   
                unitCell.expCon(9).dir=[1 4 6; 3 2 6; 2 6 6];
                unitCell.perCon=[1 4 1 9; 2 8 1 2; 4 5 1 2];
            case {'#5'}
                unitCell.Polyhedron(1)=polyhedra('rhombicuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('cube');
                unitCell.Polyhedron(3)=polyhedra('tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('tetrahedron');
                unitCell.expCon(1).dir=[3 5 1; 7 6 1; 8 14 1];
                unitCell.expCon(2).dir=[3 10 1; 2 14 1;4 22 1];
                unitCell.expCon(3).dir=[3 16 1; 2 24 1; 1 12 1];
                unitCell.perCon=[18 13; 10 7; 16 11];
            case {'#6'}
                unitCell.Polyhedron(1)=polyhedra('tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(3)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('tetrahedron');
                unitCell.expCon(1).dir=[4 1 1; 10 2 1; 2 4 1];
                unitCell.expCon(2).dir=[4 7 2; 10 8 2; 5 3 2];
                unitCell.expCon(3).dir=[3 8 3; 4 12 3; 2 11 3];
                unitCell.perCon=[1 4 1 2; 2 1 1 2; 4 2 1 2];
            case {'#7','Fig1b'}
                unitCell.Polyhedron(1)=polyhedra('cuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.expCon(1).dir=[1 3 1; 2 9 1; 3 6 1];
                unitCell.perCon=[4 6; 5 3; 1 2];
            case {'#8'}
                unitCell.Polyhedron(1)=polyhedra('truncated cube');
                unitCell.Polyhedron(2)=polyhedra('octahedron');
                unitCell.expCon(1).dir=[1 7 1; 2 9 1; 3 16 1];
                unitCell.perCon=[4 3; 6 5; 1 2];
            case {'#9','Fig6#b'}
                unitCell.Polyhedron(1)=polyhedra('rhombicuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('cuboctahedron');
                unitCell.Polyhedron(3)=polyhedra('cube');
                unitCell.Polyhedron(4)=polyhedra('cube');
                unitCell.Polyhedron(5)=polyhedra('cube');
                unitCell.expCon(1).dir=[2 14 1; 4 22 1; 1 10 1];
                unitCell.expCon(2).dir=[7 12 1; 3 24 1; 1 22 1];
                unitCell.expCon(3).dir=[7 22 1; 5 14 1; 3 21 1];
                unitCell.expCon(4).dir=[7 2 1; 3 6 1; 4 14 1];
                unitCell.perCon=[4 3;6 5;1 2];
                if strcmp(opt.template,'Fig6#b')
                    unitCell.Polyhedron(1).solidify=[1,2,3,4,5,6,8,9,11,13,16,18,19,20,21,22,24,26];
                    unitCell.Polyhedron(2).solidify=[7,9,11,12,13,14];
                    unitCell.Polyhedron(3).solidify=[3,4];
                    unitCell.Polyhedron(4).solidify=[2,5];
                    unitCell.Polyhedron(5).solidify=[2,5];
                end
            case {'#10'}
                unitCell.Polyhedron(1)=polyhedra('truncated octahedron');
                unitCell.Polyhedron(2)=polyhedra('cuboctahedron');
                unitCell.Polyhedron(3)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('truncated tetrahedron');
                unitCell.expCon(1).dir=[1 13 1; 4 19 1; 7 23 1];
                unitCell.expCon(2).dir=[6 17 1; 5 15 1; 9 23 1];
                unitCell.expCon(3).dir=[9 19 1; 10 14 1; 3 23 1];
                unitCell.perCon=[9 5 1 3;10 8 1 4; 13 8 1 3];
            case {'#11'}
                unitCell.Polyhedron(1)=polyhedra('triangular prism');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[6 4 1;5 5 1; 3 1 1];
                unitCell.perCon=[3 1 1 2;2 3 1 2;4 5 1 1];
            case {'#12','Fig1c','Fig2g'}
                unitCell.Polyhedron(1)=polyhedra('triangular prism');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[6 5 1; 3 4 1; 2 1 1];
                unitCell.expCon(2).dir=[4 5 1; 6 6 1; 3 3 1];
                unitCell.expCon(3).dir=[4 5 3; 1 6 3; 2 3 3 ];
                unitCell.perCon=[3 1 1 3;2 3 4 2;4 5 1 1];
            case {'#13'}
                unitCell.Polyhedron(1)=polyhedra('cube');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[6 8 1; 4 6 1; 1 2 1];
                unitCell.expCon(2).dir=[6 7 1; 5 5 1; 2 1 1]; 
                unitCell.perCon=[5 2 1 1;3 1 3 2;1 6 1 1];
            case {'#14','Fig6#c'}     
                unitCell.Polyhedron(1)=polyhedra('cube');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('cube');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.Polyhedron(5)=polyhedra('triangular prism');
                unitCell.Polyhedron(6)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[6 8 1; 4 6 1; 3 4 1];
                unitCell.expCon(2).dir=[5 4 2; 7 5 2; 3 2 2];                 
                unitCell.expCon(3).dir=[6 6 2; 4 5 2; 3 3 2]; 
                unitCell.expCon(4).dir=[6 5 1; 5 6 1; 3 1 1]; 
                unitCell.expCon(5).dir=[6 5 3; 5 6 3; 2 2 3]; 
                unitCell.perCon=[4 1 1 4;5 1 1 6;1 6 1 1];                
                if strcmp(opt.template,'Fig6#c')
                    unitCell.Polyhedron(1).solidify=[1 2 5 6];
                    unitCell.Polyhedron(2).solidify=[2 ];
                    unitCell.Polyhedron(4).solidify=[3 4 5];
                    unitCell.Polyhedron(5).solidify=[2];
                    unitCell.Polyhedron(6).solidify=[1];
                end
            case {'#15'}
                unitCell.Polyhedron(1)=polyhedra('triangular prism');
                unitCell.Polyhedron(2)=polyhedra('cube');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.Polyhedron(5)=polyhedra('cube');
                unitCell.Polyhedron(6)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[7 4 1; 8 5 1; 4 2 1];
                unitCell.expCon(2).dir=[6 6 2; 3 5 2; 2 1 2];
                unitCell.expCon(3).dir=[4 5 1; 6 6 1; 3 3 1];
                unitCell.expCon(4).dir=[5 6 4; 6 5 4; 1 3 4];
                unitCell.expCon(5).dir=[4 8 5; 1 7 5; 2 3 5];
                unitCell.perCon=[4 3 2 2;1 6 2 2;2 3 6 3];
            case {'#16','Fig6#d'}
                unitCell.Polyhedron(1)=polyhedra('hexagonal prism');
                unitCell.Polyhedron(2)=polyhedra('cube');
                unitCell.Polyhedron(3)=polyhedra('cube'); 
                unitCell.Polyhedron(4)=polyhedra('cube'); 
                unitCell.Polyhedron(5)=polyhedra('triangular prism'); 
                unitCell.Polyhedron(6)=polyhedra('triangular prism'); 
                unitCell.expCon(1).dir=[7 7 1; 8 8 1; 3 1 1];
                unitCell.expCon(2).dir=[7 9 1; 5 8 1; 3 3 1]; 
                unitCell.expCon(3).dir=[7 10 1; 5 9 1; 3 4 1];
                unitCell.expCon(4).dir=[4 5 4; 6 6 4; 3 2 4]; 
                unitCell.expCon(5).dir=[6 5 3; 5 6 3; 2 2 3];
                unitCell.perCon=[8 3 1 4;7 3 1 3;1 2 1 1];                
                if strcmp(opt.template,'Fig6#d')
                    unitCell.Polyhedron(1).solidify=[1 2 4 6 8];
                    unitCell.Polyhedron(2).solidify=[2 3];
                    unitCell.Polyhedron(3).solidify=[2 4];
                    unitCell.Polyhedron(4).solidify=[3 5];
                    unitCell.Polyhedron(6).solidify=[1 2 3 4 5];
                end
            case {'#17'}
                unitCell.Polyhedron(1)=polyhedra('hexagonal prism');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.Polyhedron(4)=polyhedra('triangular prism');
                unitCell.Polyhedron(5)=polyhedra('triangular prism');
                unitCell.Polyhedron(6)=polyhedra('triangular prism');
                unitCell.Polyhedron(7)=polyhedra('triangular prism');
                unitCell.Polyhedron(8)=polyhedra('triangular prism');
                unitCell.Polyhedron(9)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[4 11 1; 5 10 1; 2 4 1];
                unitCell.expCon(2).dir=[4 10 1; 5 9 1; 2 3 1];                 
                unitCell.expCon(3).dir=[4 5 2; 6 6 2; 3 3 2]; 
                unitCell.expCon(4).dir=[4 4 3; 5 6 3; 2 3 3]; 
                unitCell.expCon(5).dir=[5 4 2; 6 6 2; 3 3 2]; 
                unitCell.expCon(6).dir=[5 12 1; 6 11 1; 2 6 1];
                unitCell.expCon(7).dir=[5 5 6; 6 4 6; 3 1 6];
                unitCell.expCon(8).dir=[6 9 1; 4 8 1; 1 2 1];
                unitCell.perCon=[8 2 1 5;3 3 1 6;1 2 1 1];
            case {'#18','Fig1a','Fig2f'}
                unitCell.Polyhedron(1)=polyhedra('hexagonal prism');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[4 10 1; 6 11 1; 3 5 1];
                unitCell.expCon(2).dir=[6 7 1; 5 8 1; 2 2 1];
                unitCell.perCon=[1 8 2 1;2 4 2 1;1 2 1 1];
            case {'#19'}
                unitCell.Polyhedron(1)=polyhedra('decagonal prism');
                unitCell.Polyhedron(2)=polyhedra('triangular prism');
                unitCell.Polyhedron(3)=polyhedra('triangular prism');
                unitCell.expCon(1).dir=[1 23 1;2 24 1;4 11 1];
                unitCell.expCon(2).dir=[1 9 1; 4 21 1; 3 10 1];
                unitCell.perCon=[6 12; 10 4; 1 2];
            case {'#20'}
                unitCell.Polyhedron(1)=polyhedra('rhombicuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('octagonal prism');
                unitCell.Polyhedron(3)=polyhedra('octagonal prism');
                unitCell.Polyhedron(4)=polyhedra('octagonal prism');
                unitCell.Polyhedron(5)=polyhedra('cube');
                unitCell.Polyhedron(6)=polyhedra('cube');
                unitCell.Polyhedron(7)=polyhedra('cube');
                unitCell.Polyhedron(8)=polyhedra('truncated cube');
                unitCell.expCon(1).dir=[12 22 1; 6 14 1; 5 13 1];
                unitCell.expCon(2).dir=[6 12 1; 5 10 1; 11 22 1];
                unitCell.expCon(3).dir=[12 10 1; 4 14 1; 3 6 1];
                unitCell.expCon(4).dir=[3 4 1; 4 12 1; 1 2 1;];
                unitCell.expCon(5).dir=[3 5 1; 7 6 1; 8 14 1];
                unitCell.expCon(6).dir=[3 23 1; 1 21 1; 7 24 1];
                unitCell.expCon(7).dir=[22 9 3; 21 5 3; 13 15 3];
                unitCell.perCon=[17 8 1 3;9 4 1 2;8 8 1 2];      
            case {'#21'} 
                unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated cube');
                unitCell.Polyhedron(3)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(4)=polyhedra('truncated tetrahedron');
                unitCell.expCon(1).dir=[20 28 1; 18 26 1; 4 4 1 ];
                unitCell.expCon(2).dir=[10 48 1; 9 22 1; 5 12 1 ];
                unitCell.expCon(3).dir=[2 40 1; 12 30 1; 6 26 1];
                unitCell.perCon=[4 6; 3 5; 12 11];
            case {'#22'}
                unitCell.Polyhedron=polyhedra('cube');
                unitCell.perCon=[6 1; 4 3; 5 2];
            case {'#23'}
                unitCell.Polyhedron(1)=polyhedra('decagonal prism');
                unitCell.Polyhedron(2)=polyhedra('cube');
                unitCell.Polyhedron(3)=polyhedra('cube');
                unitCell.Polyhedron(4)=polyhedra('cube'); 
                unitCell.Polyhedron(5)=polyhedra('hexagonal prism');
                unitCell.Polyhedron(6)=polyhedra('hexagonal prism');
                unitCell.expCon(1).dir=[7 24 1; 8 13 1; 4 1 1];
                unitCell.expCon(2).dir=[7 15 1; 5 14 1; 3 3 1];
                unitCell.expCon(3).dir=[7 17 1; 5 16 1; 3 5 1];
                unitCell.expCon(4).dir=[11 14 1; 12 13 1; 5 2 1];
                unitCell.expCon(5).dir=[11 16 1; 12 15 1; 6 3 1];
                unitCell.perCon=[3 12 4 1;3 10 3 1;1 2 1 1];
            case {'#24','Fig6#a'}
                unitCell.Polyhedron(1)=polyhedra('octagonal prism');
                unitCell.Polyhedron(2)=polyhedra('cube');
                unitCell.expCon(1).dir=[7 2 1; 8 14 1; 3 1 1];
                unitCell.perCon=[4 8;6 10;1 2];                 
                if strcmp(opt.template,'Fig6#a')
                    unitCell.Polyhedron(1).solidify=[1,2,3,5,7];
                    unitCell.Polyhedron(2).solidify=[2,3,4];
                end
            case {'#25'}
                unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated octahedron');
                unitCell.Polyhedron(3)=polyhedra('cube');
                unitCell.Polyhedron(4)=polyhedra('cube');
                unitCell.Polyhedron(5)=polyhedra('cube');
                unitCell.expCon(1).dir=[1 12 1; 2 18 1; 10 48 1];
                unitCell.expCon(2).dir=[7 4 1; 3 5 1; 8 12 1];
                unitCell.expCon(3).dir=[7 22 1; 8 48 1; 3 21 1];
                unitCell.expCon(4).dir=[3 20 1; 4 43 1; 1 18 1];
                unitCell.perCon=[16 15;17 18;13 14];
            case {'#26'}
                unitCell.Polyhedron=polyhedra('hexagonal prism');
                unitCell.perCon=[7 4; 6 3; 1 2];
                opt.AZ=-34;
            case {'#27'}
                unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
                unitCell.Polyhedron(2)=polyhedra('octagonal prism');
                unitCell.Polyhedron(3)=polyhedra('octagonal prism');
                unitCell.Polyhedron(4)=polyhedra('octagonal prism');
                unitCell.expCon(1).dir=[5 43 1; 9 41 1; 7 42 1];
                unitCell.expCon(2).dir=[3 13 1; 11 5 1; 1 6 1];
                unitCell.expCon(3).dir=[5 26 1; 11 28 1; 7 20 1];
                unitCell.perCon=[24 21;26 19;23 22];
            case {'#28','Fig4a'}            
                unitCell.Polyhedron=polyhedra('truncated octahedron');
                unitCell.perCon=[14 11; 9 8; 10 7];                
                if strcmp(opt.template,'Fig4a')
                    unitCell.Polyhedron.solidify=[1,2,3,4,5,6,12,13];
                end    
            case {'#29'}
                unitCell.Polyhedron=polyhedra('rhombic dodecahedron');
                unitCell.perCon=[2 4;3 1;9 7];
            case {'#6b'}
                unitCell.Polyhedron(1)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated tetrahedron');
                unitCell.expCon(1).dir=[4 7 2; 10 8 2; 5 3 2];
                unitCell.perCon=[8 6 1 2; 7 8 1 2; 6 5 1 2];
            case {'#6b1D'}
                unitCell.Polyhedron(1)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated tetrahedron');
                unitCell.expCon(1).dir=[4 7 2; 10 8 2; 5 3 2];
                unitCell.perCon=[8 6 1 2];
            case {'#6b2D'}
                unitCell.Polyhedron(1)=polyhedra('truncated tetrahedron');
                unitCell.Polyhedron(2)=polyhedra('truncated tetrahedron');
                unitCell.expCon(1).dir=[4 7 2; 10 8 2; 5 3 2];
                unitCell.perCon=[8 6 1 2; 6 5 1 2];
            case {'#7b'}
                unitCell.Polyhedron(1)=polyhedra('cuboctahedron');
                unitCell.perCon=[4 6; 5 3; 1 2];
            case {'#7b2D'}
                unitCell.Polyhedron(1)=polyhedra('cuboctahedron');
                unitCell.perCon=[5 3; 1 2];
            case {'#7b1D'}
                unitCell.Polyhedron(1)=polyhedra('cuboctahedron');
                unitCell.perCon=[1 2];
            case {'#22_1D'}
                unitCell.Polyhedron=polyhedra('cube');
                unitCell.perCon=[6 1];
            case {'#22_2D'}
                unitCell.Polyhedron=polyhedra('cube');
                unitCell.perCon=[6 1; 5 2];
            case {'#26_1Da'}
                unitCell.Polyhedron=polyhedra('hexagonal prism');
                unitCell.perCon=[1 2];
            case {'#26_1Db'}
                unitCell.Polyhedron=polyhedra('hexagonal prism');
                unitCell.perCon=[7 4];
            case {'#26_2Da'}
                unitCell.Polyhedron=polyhedra('hexagonal prism');
                unitCell.perCon=[7 4; 6 3];
            case {'#26_2Db'}
                unitCell.Polyhedron=polyhedra('hexagonal prism');
                unitCell.perCon=[7 4; 1 2];
        end
end

% %Initialize empty array not used
% unitCell.perCon=[];
% unitCell.expCon=[];
% %Load polyhedron
% unitCell.Polyhedron(1)=polyhedra(opt.template);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE EMPTY ARRAYS THAT HAVE NOT BEEN DEFINED BEFORE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ne=1:length(unitCell.Polyhedron)
    if isfield(unitCell.Polyhedron(ne),'solidify')==0
        unitCell.Polyhedron(ne).solidify=[];
    end    
    unitCell.Polyhedron(ne).extrude=setdiff(1:length(unitCell.Polyhedron(ne).face),unitCell.Polyhedron(ne).solidify);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET OPTIONS REGARDING PERIODIC BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TURN PERIODIC BC OFF IF NO PERIODIC FACE PAIRS HAVE BEEN DEFINED
if isempty(unitCell.perCon)
    opt.periodic='off';
end
%WHEN PERIODIC FACE PAIRS ARE DEFINED, BUT OPTIONS ARE SET TO NO PERIODIC
%BC, REMOVE DEFINITION OF PERIODIC FACE PAIRS
if strcmp(opt.periodic,'off')
    for i=1:length(unitCell.Polyhedron)
        unitCell.Polyhedron(i).faceExL=ones(length(unitCell.Polyhedron(i).extrude),1);
    end
    unitCell.perCon=[];
    unitCell.perStr=[];
end
%WHEN NU SPECIFIC POLYHEDRA NUMBER IS DEFINED FOR PERIODIC FACE PAIRS, USE
%THE DEFAULT FIRST POLYHEDRON
if size(unitCell.perCon,2)==2
    unitCell.perCon=[unitCell.perCon ones(size(unitCell.perCon,1),size(unitCell.perCon,2))];
end

function [unitCell]=expandUnitCell(unitCell,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALTER POSITION OF POLYHEDRA AS DEFINED BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ne=1:(length(unitCell.Polyhedron)-1)
    centroidA=zeros(3,1); centroidB=zeros(3,1);
    for i=1:size(unitCell.expCon(ne).dir,1)
        pB(:,i)=unitCell.Polyhedron(unitCell.expCon(ne).dir(i,3)).node(unitCell.expCon(ne).dir(i,2),:)';
        pA(:,i)=unitCell.Polyhedron(ne+1).node(unitCell.expCon(ne).dir(i,1),:)';
        centroidB=centroidB+pB(:,i)/size(unitCell.expCon(ne).dir,1);
        centroidA=centroidA+pA(:,i)/size(unitCell.expCon(ne).dir,1);
    end
    H=zeros(3);
    for i=1:size(unitCell.expCon(ne).dir,1)
        H=H+(pA(:,i)-centroidA)*(pB(:,i)-centroidB)';
    end
    [U,S,V]=svd(H);
    R=V*U';
    if det(R)<0
        V(:,3)=V(:,3)*-1;
        R=V*U';
    end
    t=-R*centroidA+centroidB;
    centroidB=ones(size(unitCell.Polyhedron(ne+1).node,1),1)*centroidB';
    for i=1:size(unitCell.Polyhedron(ne+1).node,1)
        unitCell.Polyhedron(ne+1).node(i,:)=(R*unitCell.Polyhedron(ne+1).node(i,:)'+t)';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIND CONGRUENT (PERIODICALLY PLACED) FACES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIATE
faces=[];
for i=1:length(unitCell.Polyhedron)
    faces=[faces; [unitCell.Polyhedron(i).extrude' i*ones(length(unitCell.Polyhedron(i).extrude),1)]];
end
internalFacePairs=[];
externalPeriodicFacePairs=[];
unitCell.lhat=[];
unitCell.l=[];
nref=0;
rep=0;

%LOOP OVER ALL FACES ADN CHECK IF TWO FACES ARE OVERLAPPING IN THE UNIT CELL 
for i=1:size(faces,1)
    for j=(i+1):size(faces,1)
        if i~=j
            coori=unitCell.Polyhedron(faces(i,2)).node(unitCell.Polyhedron(faces(i,2)).face{faces(i,1)},:);
            coorj=unitCell.Polyhedron(faces(j,2)).node(unitCell.Polyhedron(faces(j,2)).face{faces(j,1)},:);
            if size(coori,1)==size(coorj,1)
                cooria=sum(coori)/size(coori,1);
                coorja=sum(coorj)/size(coorj,1);
                if (sum((cooria-coorja).^2)<1e-4)==1
                    st=0;
                    for k=1:size(coori,1)
                        x=coorj(:,1)-coori(k,1);
                        y=coorj(:,2)-coori(k,2); 
                        z=coorj(:,3)-coori(k,3);
                        r2=x.^2+y.^2+z.^2;
                        if sum(r2<1e-4)==0
                            st=1;
                        end
                    end
                    if st~=1
                        rep=rep+1;
                        internalFacePairs(rep,:)=[faces(i,:) faces(j,:)];
                    end
                end
            end
        end
    end
end

%LOOP OVER ALL FACES AND CHECK IF TWO FACES ARE PERIODICALLY LOCATED
if strcmp(opt.periodic,'on')
    nref=size(unitCell.perCon,1);
    %DETERMINE INITIAL LATTICE VECTORS
    for i=1:nref
        f1=unitCell.perCon(i,1);
        f2=unitCell.perCon(i,2);
        [n1, coor1]=findNormal(unitCell.Polyhedron(unitCell.perCon(i,3)),f1);
        [n2, coor2]=findNormal(unitCell.Polyhedron(unitCell.perCon(i,4)),f2);
        unitCell.lhat(i,:)=coor2-coor1;
    end

    %POSSIBLE COMBINATIONS OF LATTICE VECTORS FOR PERIODICALLY LOCATED
    %FACES
    switch nref
        case 0
            unitCell.possibleAlpha=[];
        case 1
            unitCell.possibleAlpha=1;
        case 2
            unitCell.possibleAlpha=[1 0; 0 1; 1 1; -1 1];
        case 3
            unitCell.possibleAlpha=[1 0 0; 0 1 0; 0 0 1; 1 -1 0; 1 1 0; 1 0 -1; 1 0 1; 0 -1 1; 0 1 1; 1 1 1; 1 1 -1; 1 -1 -1; -1 1 -1];
    end

    rep=0;
    for i=1:size(faces,1)
        for j=1:size(faces,1)
            if i~=j
                coori=unitCell.Polyhedron(faces(i,2)).node(unitCell.Polyhedron(faces(i,2)).face{faces(i,1)},:);
                coorj=unitCell.Polyhedron(faces(j,2)).node(unitCell.Polyhedron(faces(j,2)).face{faces(j,1)},:);          
                if size(coori,1)==size(coorj,1)
                    for l=1:size(unitCell.possibleAlpha,1)
                        coorpbc=unitCell.lhat'*unitCell.possibleAlpha(l,:)';
                        cooria=sum(coori)/size(coori,1);
                        coorja=sum(coorj)/size(coorj,1);
                        if (sum((coorja'-cooria'-coorpbc).^2)<1e-4)==1
                            st=0;
                            for k=1:size(coori,1)
                                x=coorj(:,1)-coori(k,1)-coorpbc(1); 
                                y=coorj(:,2)-coori(k,2)-coorpbc(2); 
                                z=coorj(:,3)-coori(k,3)-coorpbc(3);
                                r2=x.^2+y.^2+z.^2;
                                if sum(r2<1e-4)==0
                                    st=1;
                                end
                            end  
                            if st~=1
                                rep=rep+1;
                                externalPeriodicFacePairs(rep,:)=[faces(i,:) faces(j,:) unitCell.possibleAlpha(l,:)];
                            end
                        end
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET UP CONSTRAINTS FOR SEPARATION PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nC=size(externalPeriodicFacePairs,1);
nCi=size(internalFacePairs,1);
nP=length(unitCell.Polyhedron);

%NORMAL TO FACE PAIRS
internalFacePairsNormal=zeros(1,3);
for i=1:size(internalFacePairs)
    internalFacePairsNormal(i,:)=findNormal(unitCell.Polyhedron(internalFacePairs(i,2)),internalFacePairs(i,1));
end          

%CONSTRAINTS EXTERNAL
A=zeros(3*nC+3*nCi,nref*3+nC+3*(nP-1)+nCi);
n=zeros(3*nC+3*nCi,1);
for i=1:nC
    %Add the terms for the change in lattice vectors
    for j=1:nref
        A(3*i-2:3*i,((j-1)*3+1):(j*3))=eye(3)*externalPeriodicFacePairs(i,4+j);
    end    
    %Add the beta times normal term to the constraint matrix
    %for the periodically placed faces
    [n1 coorAve1]=findNormal(unitCell.Polyhedron(externalPeriodicFacePairs(i,2)),externalPeriodicFacePairs(i,1));
    A(3*i-2:3*i,3*nref+3*(nP-1)+i)=n1;
    %Add gamma times normal term to the constraint matrix for
    %the faces between two polyhedra in the initial template
    nf1=externalPeriodicFacePairs(i,2);
    nf2=externalPeriodicFacePairs(i,4);
    if nf1~=1
        A(3*i-2:3*i,(3*nref+3*(nf1-1)-2):(3*nref+3*(nf1-1)))=A(3*i-2:3*i,(3*nref+3*(nf1-1)-2):(3*nref+3*(nf1-1)))+eye(3);
    end
    if nf2~=1
        A(3*i-2:3*i,(3*nref+3*(nf2-1)-2):(3*nref+3*(nf2-1)))=A(3*i-2:3*i,(3*nref+3*(nf2-1)-2):(3*nref+3*(nf2-1)))-eye(3);
    end
end

%CONSTRAINTS INTERNAL FACE PAIRS
for i=1:nCi
    [n1 coorAve1]=findNormal(unitCell.Polyhedron(internalFacePairs(i,2)),internalFacePairs(i,1));
    A(3*nC+3*i-2:3*nC+3*i,3*nref+3*(nP-1)+nC+i)=n1;
    nf1=internalFacePairs(i,2);
    nf2=internalFacePairs(i,4);
    if nf1~=1
        A(3*nC+3*i-2:3*nC+3*i,(3*nref+3*(nf1-1)-2):(3*nref+3*(nf1-1)))=A(3*nC+3*i-2:3*nC+3*i,(3*nref+3*(nf1-1)-2):(3*nref+3*(nf1-1)))+eye(3);
    end
    if nf2~=1
        A(3*nC+3*i-2:3*nC+3*i,(3*nref+3*(nf2-1)-2):(3*nref+3*(nf2-1)))=A(3*nC+3*i-2:3*nC+3*i,(3*nref+3*(nf2-1)-2):(3*nref+3*(nf2-1)))-eye(3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLVE MINIMIZATION PROBLEM TO FIND EXTRUSION LENGTHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitCell.Gamma=size(A,2)-rank(A);
options=optimset('Display','off','tolfun',1e-15,'tolx',1e-15,'tolCon',1e-13,'gradObj','on');
u0=randn(size(A,2),1)+ones(size(A,2),1);
unitCell.l=[];
if size(A,1)~=0
    A=rref(A);
    [ext fval exitflag]=fmincon(@(u) ExtrudeObjective(u,nref,nP,opt),u0,[],[],A,n,[],[],[],options);
    uPol=ext(3*nref+1:3*nref+3*(nP-1));
    beta=ext(3*nref+3*(nP-1)+1:3*nref+3*(nP-1)+nC);
    gamma=ext(end-nCi+1:end);
    if exitflag<1
        fprintf('Automatic determation of the extrusion lengths\nfailed, please directly define the extrusion\nlength\n')  
    end            
    %Write extruded lentgth in correct format
    for i=1:length(unitCell.Polyhedron)
        unitCell.Polyhedron(i).faceExL=ones(1,length(unitCell.Polyhedron(i).extrude));
    end
    for i=1:size(externalPeriodicFacePairs,1)
        [val, index1]=intersect(unitCell.Polyhedron(externalPeriodicFacePairs(i,2)).extrude,externalPeriodicFacePairs(i,1));
        [val, index2]=intersect(unitCell.Polyhedron(externalPeriodicFacePairs(i,4)).extrude,externalPeriodicFacePairs(i,3));
        unitCell.Polyhedron(externalPeriodicFacePairs(i,2)).faceExL(index1)=beta(i);
        unitCell.Polyhedron(externalPeriodicFacePairs(i,4)).faceExL(index2)=beta(i);
    end
    for i=1:size(internalFacePairs,1)
        [val, index1]=intersect(unitCell.Polyhedron(internalFacePairs(i,2)).extrude,internalFacePairs(i,1));
        [val, index2]=intersect(unitCell.Polyhedron(internalFacePairs(i,4)).extrude,internalFacePairs(i,3));
        unitCell.Polyhedron(internalFacePairs(i,2)).faceExL(index1)=gamma(i);
        unitCell.Polyhedron(internalFacePairs(i,4)).faceExL(index2)=gamma(i);
    end
    %Extruded lattice vectors
    for i=1:nref
        unitCell.l(i,:)=unitCell.lhat(i,:)+2*ext(3*i-2:3*i)';
    end
    %Coordinates of polyhedra templates in extruded configuration
    for ne=1:length(unitCell.Polyhedron)
        if ne>1
            unitCell.Polyhedron(ne).nodeNew(:,1)=unitCell.Polyhedron(ne).node(:,1)+2*uPol(3*(ne-1)-2);
            unitCell.Polyhedron(ne).nodeNew(:,2)=unitCell.Polyhedron(ne).node(:,2)+2*uPol(3*(ne-1)-1);
            unitCell.Polyhedron(ne).nodeNew(:,3)=unitCell.Polyhedron(ne).node(:,3)+2*uPol(3*(ne-1));
        else
            unitCell.Polyhedron(ne).nodeNew=unitCell.Polyhedron(ne).node;
        end
    end
else
    for ne=1:length(unitCell.Polyhedron)
        unitCell.Polyhedron(ne).faceExL=ones(length(unitCell.Polyhedron(ne).extrude),1);
        unitCell.Polyhedron(ne).nodeNew=unitCell.Polyhedron(ne).node;
    end
end
unitCell.externalPeriodicFacePairs=externalPeriodicFacePairs;
unitCell.internalFacePairs=internalFacePairs;

function [n coorAve]=findNormal(unitCell,faceID)
    %Determine normal
    coor1=unitCell.node(unitCell.face{faceID}(1),:);
    coor2=unitCell.node(unitCell.face{faceID}(2),:);
    coor3=unitCell.node(unitCell.face{faceID}(3),:);
    n=cross(coor2-coor1,coor3-coor1);
    n=n/norm(n);
    %Determine face average coordinates
    coorAve=sum(unitCell.node(unitCell.face{faceID},:))/length(unitCell.face{faceID});
    
function [f df]=ExtrudeObjective(u,nref,nP,opt)
f=sum((u([(3*nref+3*(nP-1)+1):end])-opt.Lextrude).^2);
df=zeros(length(u),1);
df([(3*nref+3*(nP-1)+1):end])=2*(u([(3*nref+3*(nP-1)+1):end])-opt.Lextrude);    

function [extrudedUnitCell,unitCell]=extrudeUnitCell(unitCell,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extrudedUnitCell.edge=[];
extrudedUnitCell.edgeHinge=[];
extrudedUnitCell.node=[];
extrudedUnitCell.face=[];
extrudedUnitCell.solidify=[];
extrudedUnitCell.diagonals = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD EXTRUDED UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ne=1:length(unitCell.Polyhedron)
    nNodeSolid=size(extrudedUnitCell.node,1);
    extrudedUnitCell.node=[extrudedUnitCell.node; unitCell.Polyhedron(ne).nodeNew];
    if ~isempty(unitCell.Polyhedron(ne).solidify)
        for i=1:length(unitCell.Polyhedron(ne).solidify)
            extrudedUnitCell.solidify=[extrudedUnitCell.solidify length(extrudedUnitCell.face)+1];
            nodeNum=nNodeSolid+unitCell.Polyhedron(ne).face{unitCell.Polyhedron(ne).solidify(i)};
            extrudedUnitCell.face{end+1}=nodeNum;
            extrudedUnitCell.edgeHinge=[extrudedUnitCell.edgeHinge; nodeNum(1) nodeNum(2) nodeNum(end)];
            extrudedUnitCell.edge = [extrudedUnitCell.edge; nodeNum(1) nodeNum(2)]; 
            for j=2:length(nodeNum)-1
                extrudedUnitCell.edgeHinge=[extrudedUnitCell.edgeHinge; nodeNum(j) nodeNum(j+1) nodeNum(j-1)];
                extrudedUnitCell.edge = [extrudedUnitCell.edge; nodeNum(j) nodeNum(j+1)]; 
            end
            extrudedUnitCell.edgeHinge=[extrudedUnitCell.edgeHinge; nodeNum(end) nodeNum(1) nodeNum(end-1)];
            extrudedUnitCell.edge = [extrudedUnitCell.edge; nodeNum(end) nodeNum(1)]; 
            diagonals = nchoosek(nodeNum,2);
            edges = or(ismember(diagonals,extrudedUnitCell.edge, 'rows'),...
                ismember(fliplr(diagonals),extrudedUnitCell.edge, 'rows'));
            diagonals (edges,:) = [];
            firstdiagonal = size(extrudedUnitCell.edge,1)+1;
            extrudedUnitCell.edge = [extrudedUnitCell.edge; diagonals];
            extrudedUnitCell.diagonals = [extrudedUnitCell.diagonals firstdiagonal:size(extrudedUnitCell.edge,1)];
        end
    end
    rep=length(extrudedUnitCell.face);   
    for i=1:length(unitCell.Polyhedron(ne).extrude)        
        nNo=size(extrudedUnitCell.node,1);
        nodeNum=nNodeSolid+unitCell.Polyhedron(ne).face{unitCell.Polyhedron(ne).extrude(i)};
        nodeNumNew=nNo+1:nNo+length(nodeNum);
        unitCell.Polyhedron(ne).faceNodeExtrude{unitCell.Polyhedron(ne).extrude(i)}=nodeNumNew;
        a=extrudedUnitCell.node(nodeNum(2),:)-extrudedUnitCell.node(nodeNum(1),:);
        b=extrudedUnitCell.node(nodeNum(3),:)-extrudedUnitCell.node(nodeNum(1),:);
        alpha=acos(sum(a.*b)/(norm(a)*norm(b)));
        normal=cross(a,b)/(norm(a)*norm(b)*sin(alpha));
        extrudedUnitCell.node(nodeNumNew,:)=extrudedUnitCell.node(nodeNum,:)+...
            unitCell.Polyhedron(ne).faceExL(i)*ones(length(nodeNum),1)*normal;

        index=[1:length(nodeNum) 1];
        for i=1:length(nodeNum)   
            rep=rep+1;
            extrudedUnitCell.face{rep}=[nodeNum(index(i)) nodeNum(index(i+1)) nodeNumNew(index(i+1)) nodeNumNew(index(i))];
            extrudedUnitCell.edge([end+1],:)=[nodeNum(index(i)) nodeNum(index(i+1))];
            extrudedUnitCell.edge([end+1],:)=[nodeNum(index(i+1)) nodeNumNew(index(i+1))];
            extrudedUnitCell.edge([end+1],:)=[nodeNumNew(index(i+1)) nodeNumNew(index(i))];
            extrudedUnitCell.edge([end+1],:)=[nodeNumNew(index(i)) nodeNum(index(i))];
            extrudedUnitCell.edge([end+1],:)=[nodeNumNew(index(i)) nodeNum(index(i+1))];
            extrudedUnitCell.diagonals([end+1]) = size(extrudedUnitCell.edge,1); 
            extrudedUnitCell.edge([end+1],:)=[nodeNumNew(index(i+1)) nodeNum(index(i))];
            extrudedUnitCell.diagonals([end+1]) = size(extrudedUnitCell.edge,1);
            extrudedUnitCell.edgeHinge([end+1],:)=[nodeNum(index(i)) nodeNum(index(i+1)) nodeNumNew(index(i+1)) nodeNumNew(index(i))];
            extrudedUnitCell.edgeHinge([end+1],:)=[nodeNum(index(i+1)) nodeNumNew(index(i+1)) nodeNumNew(index(i)) nodeNum(index(i))];
            extrudedUnitCell.edgeHinge([end+1],:)=[nodeNumNew(index(i)) nodeNum(index(i)) nodeNum(index(i+1)) nodeNumNew(index(i+1))]; 
        end
    end
    nodes(ne) = size(extrudedUnitCell.node,1);
end
extrudedUnitCell.edgeHingeSorted=extrudedUnitCell.edgeHinge;
extrudedUnitCell.edgeHingeSorted(:,1:2)=sort(extrudedUnitCell.edgeHinge(:,1:2),2);

[a,b,c] = unique(extrudedUnitCell.edgeHingeSorted(:,1:2),'rows');
A=[a accumarray(c,1)];
sharedEdge=A(A(:,3)==2,1:2);
nodeHingeEx=zeros(size(sharedEdge,1),4);
% nodeHingeEx=zeros(size(sharedEdge,1)*2,4);

for i=1:size(sharedEdge,1)
    refNode=extrudedUnitCell.edgeHingeSorted(((extrudedUnitCell.edgeHingeSorted(:,1)==sharedEdge(i,1)).*...
        (extrudedUnitCell.edgeHingeSorted(:,2)==sharedEdge(i,2)))==1,3)';
%     refNode=extrudedUnitCell.edgeHingeSorted(((extrudedUnitCell.edgeHingeSorted(:,1)==sharedEdge(i,1)).*...
%         (extrudedUnitCell.edgeHingeSorted(:,2)==sharedEdge(i,2)))==1,3:4)';
    %Determine order
    refNode2=extrudedUnitCell.edgeHinge(((extrudedUnitCell.edgeHinge(:,1)==sharedEdge(i,1)).*...
        (extrudedUnitCell.edgeHinge(:,2)==sharedEdge(i,2)))==1,3)';

    if length(refNode2)>1
        'error'
    end
    index=1:2;
    index1=index(refNode(1,:)==refNode2(1));
    index2=setdiff(index,index1);
    nodeHingeEx(i,:)=[sharedEdge(i,:) refNode(1,[index1,index2])];
%     nodeHingeEx(2*i-1,:)=[sharedEdge(i,:) refNode(1,[index1,index2])];
%     nodeHingeEx(2*i,:)=[sharedEdge(i,:) refNode(2,[index1,index2])];
end
extrudedUnitCell.nodeHingeEx=nodeHingeEx;

%Determine initial edge length
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    extrudedUnitCell.edgeL(i)=sqrt(dx*dx');
end

%Determine initial angles
extrudedUnitCell.theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [~,extrudedUnitCell.theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end

%Reference Nodes
if strcmp(opt.periodic,'on')
    extrudedUnitCell.ref=[];
    nref=size(unitCell.l,1);
    for i=1:nref
        extrudedUnitCell.node(end+1,:)=0.01*randn(3,1); 
        extrudedUnitCell.ref(i)=size(extrudedUnitCell.node,1);
    end
end

%Internal Hinges
%added by Agustin Iniguez
prevPolNodes = 0;
i = 1;
for ne=1:length(unitCell.Polyhedron)
    internalHinges = length(unitCell.Polyhedron(ne).node)+prevPolNodes;
    for hinge = 1: length(extrudedUnitCell.nodeHingeEx)
        if (extrudedUnitCell.nodeHingeEx(hinge,1) <= internalHinges && ...
                extrudedUnitCell.nodeHingeEx(hinge,2) <= internalHinges && ...
                extrudedUnitCell.nodeHingeEx(hinge,1) > prevPolNodes &&...
                extrudedUnitCell.nodeHingeEx(hinge,2) > prevPolNodes)
            extrudedUnitCell.innerHinges(i) = hinge;
            i = i+1;
        end
    end
    prevPolNodes = nodes(ne);
end

%Max stretch possible
if strcmp(opt.constrEdge,'off')
    if ~isnan(opt.maxStretch)
        extrudedUnitCell.maxStretch = sum(extrudedUnitCell.edgeL*opt.maxStretch.^2);
    else
        extrudedUnitCell.maxStretch = inf;
    end
else
    extrudedUnitCell.maxStretch = 0;
end

%Max hinge folding possible
if ~isnan(opt.constAnglePerc)
    extrudedUnitCell.maxHingeFold = sum(max(abs(pi*opt.constAnglePerc+extrudedUnitCell.theta),...
        abs(pi*opt.constAnglePerc-extrudedUnitCell.theta)).^2);
else 
    extrudedUnitCell.maxHingeFold = inf;
end


