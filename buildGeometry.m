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
    case {'individual'}
        %Initialize empty array not used
        unitCell.perCon=[];
        unitCell.expCon=[];
        %Load polyhedron
        unitCell.Polyhedron(1)=polyhedra(opt.template);
    case {'material'}
        %LOAD ONE OF THE PREDEFINED UNIFORM SPACE-FILLING TESSELATIONS
        unitCell = createMaterial(opt.template);
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
%WHEN NO SPECIFIC POLYHEDRA NUMBER IS DEFINED FOR PERIODIC FACE PAIRS, USE
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

for i=1:size(sharedEdge,1)
    refNode=extrudedUnitCell.edgeHingeSorted(((extrudedUnitCell.edgeHingeSorted(:,1)==sharedEdge(i,1)).*...
        (extrudedUnitCell.edgeHingeSorted(:,2)==sharedEdge(i,2)))==1,3)';
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
end
extrudedUnitCell.nodeHingeEx=nodeHingeEx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine initial edge length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    extrudedUnitCell.edgeL(i)=sqrt(dx*dx');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine initial angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extrudedUnitCell.theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [~,extrudedUnitCell.theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Reference Nodes for lattice vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(opt.periodic,'on')
    extrudedUnitCell.ref=[];
    nref=size(unitCell.l,1);
    for i=1:nref
        extrudedUnitCell.node(end+1,:)=0.01*randn(3,1); 
        extrudedUnitCell.ref(i)=size(extrudedUnitCell.node,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Internal Hinges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Max stretch possible calculation for max energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(opt.constrEdge,'off')
    if ~isnan(opt.maxStretch)
        extrudedUnitCell.maxStretch = sum(extrudedUnitCell.edgeL*opt.maxStretch.^2);
    else
        extrudedUnitCell.maxStretch = inf;
    end
else
    extrudedUnitCell.maxStretch = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Max hinge folding possible for max energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnan(opt.constAnglePerc)
    extrudedUnitCell.maxHingeFold = sum(max(abs(pi*opt.constAnglePerc+extrudedUnitCell.theta),...
        abs(pi*opt.constAnglePerc-extrudedUnitCell.theta)).^2);
else 
    extrudedUnitCell.maxHingeFold = inf;
end


