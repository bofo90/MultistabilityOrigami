function unitCell = createMaterial(name)

switch name
    case {'rhombicuboctahedron_mat1'}
        unitCell.Polyhedron(1)=polyhedra('rhombicuboctahedron');
        unitCell.perCon=[18 13; 10 7; 16 11];
    case {'truncated tetrahedron_mat'}
        unitCell.Polyhedron(1)=polyhedra('truncated tetrahedron');
        unitCell.Polyhedron(2)=polyhedra('truncated tetrahedron');
        unitCell.expCon(1).dir=[4 7 2; 10 8 2; 5 3 2];
        unitCell.perCon=[8 6 1 2; 7 8 1 2; 6 5 1 2];
    case {'cuboctahedron_mat'}
        unitCell.Polyhedron(1)=polyhedra('cuboctahedron');
        unitCell.perCon=[4 6; 5 3; 1 2];
    case {'truncated cube_mat'}
        unitCell.Polyhedron(1)=polyhedra('truncated cube');
        unitCell.perCon=[4 3; 6 5; 1 2];
    case {'rhombicuboctahedron_mat2'}
        unitCell.Polyhedron(1)=polyhedra('rhombicuboctahedron');
        unitCell.perCon=[4 3;6 5;1 2];
    case {'triangular prism_mat1'}
        unitCell.Polyhedron(1)=polyhedra('triangular prism');
        unitCell.Polyhedron(2)=polyhedra('triangular prism');
        unitCell.expCon(1).dir=[6 4 1;5 5 1; 3 1 1];
        unitCell.perCon=[3 1 1 2;2 3 1 2;4 5 1 1];
    case {'triangular prism_mat2'}
        unitCell.Polyhedron(1)=polyhedra('triangular prism');
        unitCell.Polyhedron(2)=polyhedra('triangular prism');
        unitCell.Polyhedron(3)=polyhedra('triangular prism');
        unitCell.Polyhedron(4)=polyhedra('triangular prism');
        unitCell.expCon(1).dir=[6 5 1; 3 4 1; 2 1 1];
        unitCell.expCon(2).dir=[4 5 1; 6 6 1; 3 3 1];
        unitCell.expCon(3).dir=[4 5 3; 1 6 3; 2 3 3 ];
        unitCell.perCon=[3 1 1 3;2 3 4 2;4 5 1 1];
    case {'dodecagonal prism'}
        unitCell.Polyhedron(1)=polyhedra('dodecagonal prism');
        unitCell.perCon=[6 12; 10 4; 1 2];
    case {'truncated cuboctahedron_mat1'} 
        unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
        unitCell.perCon=[4 6; 3 5; 12 11];
    case {'cube_mat'}
        unitCell.Polyhedron=polyhedra('cube');
        unitCell.perCon=[6 1; 4 3; 5 2];
    case {'octagonal prism_mat',}
        unitCell.Polyhedron(1)=polyhedra('octagonal prism');
        unitCell.perCon=[4 8;6 10;1 2]; 
    case {'truncated cuboctahedron_mat2'}
        unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
        unitCell.perCon=[16 15;17 18;13 14];
    case {'hexagonal prism_mat'}
        unitCell.Polyhedron=polyhedra('hexagonal prism');
        unitCell.perCon=[7 4; 6 3; 1 2];
    case {'truncated cuboctahedron_mat3'}
        unitCell.Polyhedron(1)=polyhedra('truncated cuboctahedron');
        unitCell.perCon=[24 21;26 19;23 22];
    case {'truncated octahedron_mat'}            
        unitCell.Polyhedron=polyhedra('truncated octahedron');
        unitCell.perCon=[14 11; 9 8; 10 7];                
end