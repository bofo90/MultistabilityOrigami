function findDeformation(unitCell,extrudedUnitCell,opt)

%Show details geometries (if requested)
if strcmp(opt.analysis,'result')
    fprintf('Maximum stretching %1.2f.\n', opt.maxStretch);
    fprintf('kH %f\tkTA %f\tkE %f\tkF %f\n', opt.Khinge, opt.KtargetAngle, opt.Kedge, opt.Kface);
    switch opt.readHingeFile
        case 'off'
            metadataFile(opt, unitCell, extrudedUnitCell);
            nonlinearFolding(unitCell,extrudedUnitCell,opt,opt.angleConstrFinal(1).val);
%             
        case 'on'
            opt.angleConstrFinal = [];
            fileHinges = strcat(pwd, '/Results/hingeList_reduced/', opt.template, '.csv');
            if ~exist(fileHinges, 'file')
                fprintf('Hinge-selection file does not exist.\n');
            else
                hingeList = dlmread(fileHinges);
                metadataFile(opt, unitCell, extrudedUnitCell);
                maxHinge = opt.maxHinges;
                minHinge = opt.minHinges;
                numHinges = size(hingeList, 1);
                anglFold = -(pi-pi*(opt.constAnglePerc-0.005));
                parfor i = 1:numHinges
                    row = hingeList(i, :);
                    hinges = row(0~=row);
                    if length(hinges) <= maxHinge && length(hinges) >= minHinge
                        angles = [hinges(:), anglFold * ones(length(hinges), 1)];
                        fprintf('Hinge selection number %d/%d.\n', i, numHinges);
                        nonlinearFolding(unitCell,extrudedUnitCell,opt, angles);
                    end
                end
            end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NON-LINEAR ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonlinearFolding(unitCell,extrudedUnitCell,opt,angtemp)

%INITIALIZE LINEAR CONSTRAINTS
[Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt);

%Save some variables
opt.angleConstrFinal(1).val = angtemp;
extrudedUnitCell.angleConstr=[];
result = [];
E=[];
exfl= [];
u0=zeros(3*size(extrudedUnitCell.node,1),1);
theta0=extrudedUnitCell.theta;

%Create file for saving the results
extraName = sprintf('/kh%2.3f_kta%2.3f_ke%2.3f_kf%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge, opt.Kface);
folderName = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat', opt.saveFile, extraName);
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
    
%%%%%% Folding part %%%%%%
%Run the Folding of the structure
initialiseGlobalx(u0, theta0);
angles = theta0(opt.angleConstrFinal(1).val(:,1)) + linspace(0,1,opt.steps+1).*(opt.angleConstrFinal(1).val(:,2) - theta0(opt.angleConstrFinal(1).val(:,1)));
V(:,1)=u0;
[E.E(1,1),~,E.Eedge(1,1),E.Ediag(1,1),E.Eface(1,1),E.Ehinge(1,1),E.EtargetAngle(1,1), ~]=Energy(u0,extrudedUnitCell,opt);
exfl(1,1) = 1;
exfl(2,1) = 1;
Etemp = [];
for anglestep = 2:opt.steps+1
    exfltemp = [];
    opt.angleConstrFinal(1).val(:,2) = angles(:,anglestep);
    [Vtemp, exfltemp, output, ~] = FoldStructure(u0, Etemp, exfltemp, extrudedUnitCell, opt, 1, Aeq, Beq);
    u0 = Vtemp(:,2);
    if exfltemp(2,1) ~= 1
        exfl(2,1) = exfltemp(2,1);
    end
end
V(:,2)=u0;
extrudedUnitCell.angleConstr=opt.angleConstrFinal(1).val;
[E.E(2,1),~,E.Eedge(2,1),E.Ediag(2,1),E.Eface(2,1),E.Ehinge(2,1),E.EtargetAngle(2,1), ~]=Energy(u0,extrudedUnitCell,opt);
extrudedUnitCell.angleConstr=[];
[result, theta1,u1] = SaveResultPos(result, opt, V, output, 1);

%%%%%% Releasing part %%%%%%
%change algorithm for releasing
opt.options.Algorithm = opt.relAlgor;
opt.angleConstrFinal(2).val = [];
opt.KtargetAngle = 0;

initialiseGlobalx(u1, theta1);
[V, exfl, output, E] = FoldStructure(u1, E, exfl, extrudedUnitCell, opt, 2, Aeq, Beq);
[result, ~,~] = SaveResultPos(result, opt, V, output, 2);

%Return to original options
opt.options.Algorithm = opt.folAlgor;
result = SaveResultEnergy(result, E, exfl, opt);

%Save the result in a file
fileName = strcat(folderName,'/',mat2str(opt.angleConstrFinal(1).val(:,1)'),'op.mat');
save(fileName, 'result');


%Clear variables for next fold
clearvars result E exfl output;
fclose('all');

function [V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, iter, Aeq, Beq)

V = [];
V(:,1)=u0;
[E.E(1,iter),~,E.Eedge(1,iter),E.Ediag(1,iter),E.Eface(1,iter),E.Ehinge(1,iter),E.EtargetAngle(1,iter), ~]=Energy(u0,extrudedUnitCell,opt);
exfl(1,iter) = 1;


%Run the Folding of the structure
if isempty(opt.angleConstrFinal(iter).val)
    fprintf('Angle contrain: None\n');
else
    fprintf(['Angle contrain:', mat2str(opt.angleConstrFinal(iter).val(:,1)') ,'\n']);
end
extrudedUnitCell.angleConstr=opt.angleConstrFinal(iter).val;
% fprintf('Folding:\t');
% t1 = toc;
%Determine new equilibrium
[V(:,2),~,exfl(2,iter),output]=fmincon(@(u) Energy(u,extrudedUnitCell,opt),u0,[],[],Aeq,Beq,[],[],@(u) nonlinearConstr(u,extrudedUnitCell,opt),opt.options);
u0 = V(:,2);
%Determine energy of that equilibrium
[E.E(2,iter),~,E.Eedge(2,iter),E.Ediag(2,iter),E.Eface(2,iter),E.Ehinge(2,iter),E.EtargetAngle(2,iter), ~]=Energy(u0,extrudedUnitCell,opt);
% t2 = toc;
% fprintf('time %1.2f, exitflag %d\n',t2-t1,exfl(2,iter));


function [result, lastAngle, lastPosition] = SaveResultPos(result, opt, Positions, minimizationOuput, state)

%Get Angles from global variable
angles = getGlobalAngles;     
switch opt.gethistory
    case 'on'
        Positions = getGlobalAllx;
    case 'off'
        angles = [angles(:,1) angles(:,end)];
end

%Arrange results in a better way to save them
result.deform(state).V=[Positions(1:3:end,end) Positions(2:3:end,end) Positions(3:3:end,end)];
result.deform(state).Ve=Positions(:,end);
result.deform(state).theta = angles(:,end);
result.deform(state).output = minimizationOuput;

if ~isfield(result.deform(state), 'interV')
    result.deform(state).interV = [];
end

for j=1:size(Positions,2)
    result.deform(state).interV(end+1).V=[Positions(1:3:end,j) Positions(2:3:end,j) Positions(3:3:end,j)];
    result.deform(state).interV(end).Ve=Positions(:,j);
    result.deform(state).interV(end).theta = angles(:,j);
end

lastAngle = angles(:,end);
lastPosition = Positions(:,end);

function result = SaveResultEnergy(result, E, exfl, opt)

result.E=E.E;
result.Eedge=E.Eedge;
result.Ediag=E.Ediag;
result.Eface=E.Eface;
result.Ehinge=E.Ehinge;
result.EtargetAngle=E.EtargetAngle;    
result.exfl = exfl;
result.numMode=length(result.deform);
result.anglConstr = opt.angleConstrFinal(1).val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENERGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, dE,Eedge,Ediag,Eface,Ehinge,EtargetAngle,theta]=Energy(u,extrudedUnitCell,opt)

E=0; dE=zeros(3*size(extrudedUnitCell.node,1),1);
Eedge=0;
Ediag = 0;
Eface=0;
Ehinge=0;
EtargetAngle=0;

%APPLY DEFORMATION NODES
extrudedUnitCellPrev = extrudedUnitCell;
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%ENERGY ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'off')
    [dEdge, Jedge]=getEdge(extrudedUnitCell);
    %shearing energy
    Ediag=1/2*opt.Kdiag*sum(dEdge(extrudedUnitCell.diagonals).^2);
    dE=dE+opt.Kdiag*Jedge(extrudedUnitCell.diagonals,:)'*dEdge(extrudedUnitCell.diagonals);
    %streching energy
    notdiagonal = 1:size(extrudedUnitCell.edge,1);
    notdiagonal(extrudedUnitCell.diagonals) = [];
    Eedge=1/2*opt.Kedge*sum(dEdge(notdiagonal).^2);
    dE=dE+opt.Kedge*Jedge(notdiagonal,:)'*dEdge(notdiagonal);
end

%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'off')
    [dFace, Jface]=getFace(extrudedUnitCell);
    Eface=opt.Kface/2*sum(dFace.^2);
    dE=dE+opt.Kface*Jface'*dFace;
end

%ENERGY ASSOCIATED TO HINGE BENDING
[theta, Jhinge]=getHinge(extrudedUnitCell, extrudedUnitCellPrev);
Ehinge=1/2*opt.Khinge*sum((theta-extrudedUnitCell.theta).^2);
dE=dE+opt.Khinge*(Jhinge'*(theta-extrudedUnitCell.theta));

%ENERGY ASSOCIATED TO TARGET HINGE ANGLES
if size(extrudedUnitCell.angleConstr,1)==0
    dtheta=[];
else
    dtheta=theta(extrudedUnitCell.angleConstr(:,1))-extrudedUnitCell.angleConstr(:,2);
    dE=dE+opt.KtargetAngle*Jhinge(extrudedUnitCell.angleConstr(:,1),:)'*dtheta;
end
EtargetAngle=1/2*opt.KtargetAngle*sum(dtheta.^2);

%TOTAL ENERGY
E=Eedge+Ediag+Eface+Ehinge+EtargetAngle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NONLINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,Ceq,DC,DCeq]=nonlinearConstr(u,extrudedUnitCell,opt)

C1=[]; C2=[];
DC1=[]; DC2=[];
Ceq1=[]; Ceq2=[]; Ceq3=[];
DCeq1=[]; DCeq2=[]; DCeq3=[];

%APPLY DEFORMATION NODES
extrudedUnitCellPrev = extrudedUnitCell;
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%%%%%%%%%%%%%%%%%%%%%%
%INEQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM AND MINIMUM ANGLES
if ~isnan(opt.constAnglePerc)
    [angles, Dangles]=getHinge(extrudedUnitCell, extrudedUnitCellPrev);
    C1 = [-angles-pi*opt.constAnglePerc; angles-pi*opt.constAnglePerc];
    DC1 = [-Dangles; Dangles];
end
%MAXIMUM AND MINIMUM EDGE STRECHING
if strcmp(opt.constrEdge,'off') && ~isnan(opt.maxStretch)
    [normStrech, DnormStrech]=getEdgeNorm(extrudedUnitCell);
    C2 = [-normStrech-opt.maxStretch; normStrech-opt.maxStretch];
    DC2 = [-DnormStrech; DnormStrech];
end
C = [C1; C2];
DC = [DC1; DC2]';


%%%%%%%%%%%%%%%%%%%%%%
%EQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINT ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'on')
    [Ceq1, DCeq1]=getEdge(extrudedUnitCell);
end
%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'on')
    [Ceq2, DCeq2]=getFace(extrudedUnitCell);
end
%FINAL CONSTRAINTS
Ceq=[Ceq1; Ceq2; Ceq3];
DCeq=[DCeq1; DCeq2; DCeq3]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt)

%FIX NODE CONSTRAINTS
%IMPROVE FOLLOWING - AUTOMATIC DEPENDING ON NODES OF FACE 1
nodeFix=extrudedUnitCell.face{1};
e1=extrudedUnitCell.node(nodeFix(2),:)-extrudedUnitCell.node(nodeFix(1),:);
e2=extrudedUnitCell.node(nodeFix(3),:)-extrudedUnitCell.node(nodeFix(2),:);
e1=e1/norm(e1);
e2=e2/norm(e2);
e3=cross(e2,e1);
e3=e3/norm(e3);

Aeq=zeros(6,3*size(extrudedUnitCell.node,1));
Beq=zeros(6,1);
Aeq(1,3*nodeFix(2)-2)=1;
Aeq(2,3*nodeFix(2)-1)=1;
Aeq(3,3*nodeFix(2))=1;
Aeq(4,3*nodeFix(1)-2:3*nodeFix(1))=e3;
Aeq(5,3*nodeFix(1)-2:3*nodeFix(1))=e2;
Aeq(6,3*nodeFix(3)-2:3*nodeFix(3))=e3;

%MERGE NODES AT INITIALLY SAME LOCATION
rep=size(Aeq,1);
for i=1:size(unitCell.internalFacePairs,1)   
    for j=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)})
        index1=unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)}(j);
        for k=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)})
            index2=unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)}(k);
            if norm(extrudedUnitCell.node(index2,:)'-extrudedUnitCell.node(index1,:)')<opt.Lextrude/1e6
                rep=rep+1;
                Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                Aeq(3*rep-2:3*rep,3*index1-2:3*index1)=[1 0 0; 0 1 0; 0 0 1];
                Aeq(3*rep-2:3*rep,3*index2-2:3*index2)=[-1 0 0; 0 -1 0; 0 0 -1];
                Beq(3*rep-2:3*rep)=0;
            end
        end
     end
end

%PERIODIC NODAL CONSTRAINTS
if strcmp(opt.periodic,'on')
    nref=length(extrudedUnitCell.ref);
    rep=size(Aeq,1);
    for i=1:size(unitCell.possibleAlpha,1)
        for j=1:size(extrudedUnitCell.node,1)-nref
            coor1=extrudedUnitCell.node(j,:)';
            for k=1:size(extrudedUnitCell.node,1)-nref
                coor2=extrudedUnitCell.node(k,:)';
                if norm(coor2-coor1-unitCell.l'*unitCell.possibleAlpha(i,:)')<1e-6
                    rep=rep+1;
                    %sprintf('%d, node 1 = %d, node 2 =%d',[rep,j,k])
                    Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                    Aeq(3*rep-2:3*rep,3*j-2:3*j)=[1 0 0; 0 1 0; 0 0 1];
                    Aeq(3*rep-2:3*rep,3*k-2:3*k)=[-1 0 0; 0 -1 0; 0 0 -1];
                    for l=1:nref
                        Aeq(3*rep-2:3*rep,3*extrudedUnitCell.ref(l)-2:3*extrudedUnitCell.ref(l))=unitCell.possibleAlpha(i,l)*[-1 0 0; 0 -1 0; 0 0 -1];
                    end
                    Beq(3*rep-2:3*rep)=0;
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDGE LENGTH AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dEdge, Jedge]=getEdge(extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
Jedge=zeros(length(extrudedUnitCell.edge),size(extrudedUnitCell.node,1)*3);   
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=L-extrudedUnitCell.edgeL(i);            
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))=-dx/L;
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))=dx/L;
end


function [dEdge, Jedge]=getEdgeNorm(extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
Jedge=zeros(length(extrudedUnitCell.edge),size(extrudedUnitCell.node,1)*3);   
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=(L-extrudedUnitCell.edgeL(i))/extrudedUnitCell.edgeL(i);            
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))=-dx/L/extrudedUnitCell.edgeL(i);
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))=dx/L/extrudedUnitCell.edgeL(i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FACE OUT OF PLANE DEFORMATION AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dFace, Jface]=getFace(extrudedUnitCell)
rep=0;
tnf=0;
for i=1:length(extrudedUnitCell.face)
    tnf=tnf+length(extrudedUnitCell.face{i})-3;
end
dFace=zeros(tnf,1);
Jface=zeros(tnf,size(extrudedUnitCell.node,1)*3);
for i=1:length(extrudedUnitCell.face)
    coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(2),:);
    coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(3),:);
    a=cross(coor2-coor1,coor3-coor1);
    for j=1:length(extrudedUnitCell.face{i})-3
        rep=rep+1;
        coor4=extrudedUnitCell.node(extrudedUnitCell.face{i}(3+j),:);
        Jface(rep,3*extrudedUnitCell.face{i}(1)-2:3*extrudedUnitCell.face{i}(1))=cross((coor3-coor2),(coor3-coor4));
        Jface(rep,3*extrudedUnitCell.face{i}(2)-2:3*extrudedUnitCell.face{i}(2))=cross((coor3-coor1),(coor4-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3)-2:3*extrudedUnitCell.face{i}(3))=cross((coor4-coor1),(coor2-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3+j)-2:3*extrudedUnitCell.face{i}(3+j))=cross((coor2-coor1),(coor3-coor1));
        dFace(rep)=(coor4-coor1)*a';
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HINGE ANGLE AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, Jhinge]=getHinge(extrudedUnitCell, extrudedUnitCellPrev)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
thetaPrev = getGlobalx(extrudedUnitCellPrev);
Jhinge=zeros(size(extrudedUnitCell.nodeHingeEx,1),size(extrudedUnitCell.node,1)*3);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [Jhinge(i,index),theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
    if abs(thetaPrev(i)-theta(i)) >= 0.50*2*pi
        theta(i) = theta(i)+sign(thetaPrev(i))*2*pi;
%         fprintf('angle %d energy\n', i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions to extract information from Global Variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initialiseGlobalx(u0, theta)
global poop
poop.x = [];
poop.theta = [];
poop.x = [poop.x u0];
poop.theta = [poop.theta theta];
poop.flag = 0;

function theta = getGlobalx(extrudedUnitCell)
global poop
if poop.flag
    theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
    extrudedUnitCell.node=extrudedUnitCell.node+[poop.x(1:3:end,end) poop.x(2:3:end,end) poop.x(3:3:end,end)];
    for i=1:size(extrudedUnitCell.nodeHingeEx,1)
        [~,theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
        if abs(poop.theta(i,end)-theta(i)) >= 0.5*2*pi
            theta(i) = theta(i)+sign(poop.theta(i,end))*2*pi;
%             fprintf('angle %d global\n', i);
        end
    end
    poop.theta = [poop.theta theta];
    %turn the flag off after analising the angles just after an iteration,
    %so it doesnt analyse it every time you get the energy
    poop.flag = 0;
else
    theta = poop.theta(:,end);
end

function r = getGlobalAllx
global poop
r = poop.x;

function theta = getGlobalAngles
global poop
theta = poop.theta;

