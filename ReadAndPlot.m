function ReadAndPlot(unitCell, extrudedUnitCell, opt)

switch opt.analysis
    case 'info'
        result = [];
        outputResults(unitCell,extrudedUnitCell,result,opt);
    case {'result', 'savedata', 'plot'}
        %get file of results
        extraName = sprintf('/kh%2.3f_kta%2.3f_ke%2.3f_kf%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge, opt.Kface);
        folderResults = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat', opt.saveFile, extraName);
        if ~exist(folderResults, 'dir')
            fprintf(['No folder with results:',folderResults,'\n']);
        else
            
            %create folder of data with files
            if strcmp(opt.analysis,'savedata')
                folderEnergy = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/energy', opt.saveFile, extraName);
                [fMassDist, fHinge, fEnergy, fAngles] = makeFileswHeaders(folderEnergy, folderResults);
                
            end
            
            allFiles = dir(folderResults);
            directories = 0;
            succesfullFiles = 0;
            
            for ct = 1:length(allFiles)
                if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:end-4), 'metadata')
                    % skip all directories and metadata file
                    directories = directories+1;
                    continue;
                end

                % parse the file name to get back hinge set
                resfilename = allFiles(ct).name;
                [hingeSet, ~] = getHingeSet(resfilename);
                if strcmp(opt.readHingeFile,'off')
                    if ~isequal(hingeSet, opt.angleConstrFinal(end).val(:,1))
                        continue;
%                     elseif ~strcmp(resfilename(1:end-4), '[8 3]_Ang1_18_Angl2_27')
%                         continue;
                    end
                end
                extrudedUnitCell.angleConstr = [hingeSet(:), -pi*opt.constAnglePerc*ones(length(hingeSet), 1)];
                % load results from file
                load(strcat(folderResults,'/', allFiles(ct).name), 'result');
                succesfullFiles = succesfullFiles + 1;
                fprintf('Plot of Hinges number %d/%d\n', succesfullFiles, length(allFiles)-directories);
                
                if strcmp(opt.analysis, 'savedata')
                    [CM, Radios, Stdev, EhineInt, SumIntAngles, SumExtAngles, maxStrech, minStrech] =...
                                                getData(extrudedUnitCell, opt, result);
                    Energies = [ones(size(result.E,2),1)*(ct-directories), result.Eedge(2,:)',...
                        result.Ediag(2,:)', result.Eface(2,:)', result.Ehinge(2,:)',...
                        result.EtargetAngle(2,:)', EhineInt(2,:)', result.exfl(2,:)'];
                    PosStad = [ones(size(result.E,2),1,1)*(ct-directories),...
                        reshape(CM(2,:)',[2,3]),Radios(2,:)', Stdev(2,:)',...
                        maxStrech(2,:)', minStrech(2,:)', SumIntAngles(2,:)', SumExtAngles(2,:)'];
                    Hinges = [num2str(ct-directories),',',mat2str(hingeSet'),',',...
                        mat2str(result.anglConstr(1,2),5)];
                    AllAngles = [extrudedUnitCell.theta];
                    for iter = 1:size(result.deform,2)
                        AllAngles = [AllAngles result.deform(iter).theta];
                    end
                    AllAngles = [ones(size(AllAngles,2),1)*(ct-directories) AllAngles'];
                    dlmwrite(fMassDist, PosStad, 'delimiter', ',', '-append','precision',7);
                    dlmwrite(fHinge, Hinges, 'delimiter', '', '-append');
                    dlmwrite(fEnergy, Energies, 'delimiter', ',', '-append','precision',7);
                    dlmwrite(fAngles, AllAngles, 'delimiter', ',', '-append','precision',7);
                end
                
                if strcmp(opt.createFig, 'on') || strcmp(opt.analysis, 'plot')
                    
                    nameFolderPlot=[pwd,'/Results/',opt.template,'/',opt.relAlgor,'/images',...
                        opt.saveFile,extraName];
                    nameFilePlot = ['/',resfilename(1:end-4),'_AnglEv.png'];
                    
                    if ~exist(nameFolderPlot, 'dir')
                        mkdir(nameFolderPlot);
                    end
                    allangles = [];
                    for iter = 1:size(result.deform,2)
                        allangles = [allangles result.deform(iter).interV(:).theta];
                    end
                    plot(allangles')
                    x = 0;
                    for iter = 1:(size(result.deform,2)-1)
                        x = x + size(result.deform(iter).interV,2)+0.5;
                        line([x x],[-1.1*pi 1.1*pi], 'Color', [0 0 0])
                    end
                    saveas(gcf, [nameFolderPlot, nameFilePlot]);
                    savefig([nameFolderPlot,'/',resfilename(1:end-4),'_AnglEv.fig'])
                    close 'all';                    
                    
                    outputResults(unitCell,extrudedUnitCell,result,opt,resfilename(1:end-4));
                end
                close all;
            end
        end
        fclose('all');
        
end

function [fileMassDist, fileHinge, fileEnergy, fileAngles] = makeFileswHeaders(folderEnergy, folderResults)

if ~exist(folderEnergy, 'dir')
    mkdir(folderEnergy);
end

fileEnergy = strcat(folderEnergy, '/','EnergyData.csv');
if exist(fileEnergy, 'file')
    delete(fileEnergy) % always start with new file
end
headersEnergy = {'Hinge Number'; 'EdgeEnergy';'DiagonalEnergy'; 'FaceEnergy'; ...
    'HingeEnergy'; 'TargetAngleEnergy'; 'InternalHingeEnergy'; 'Flags'};
writeHeader(fileEnergy, headersEnergy);

fileHinge = strcat(folderEnergy, '/','Hinges.csv');
if exist(fileHinge, 'file')
    delete(fileHinge) % always start with new file
end
headersHinge = {'HingeNumber'; 'ActuatedHinges'; 'Theta1'};
writeHeader(fileHinge, headersHinge);

fileMassDist = strcat(folderEnergy, '/','PosStad.csv');
if exist(fileMassDist, 'file')
    delete(fileMassDist) % always start with new file
end
headersMassDist = {'Hinge Number';'CenterMassXFol';'CenterMassXRel';'CenterMassYFol';'CenterMassYRel';...
    'CenterMassZFol';'CenterMassZRel'; 'MeanDistanceCMFol';'MeanDistanceCMRel'; 'StdDevDistanceCMFol';...
    'StdDevDistanceCMRel';'MaxEdgeStrechFol';'MaxEdgeStrechRel';'MinEdgeStrechFol';'MinEdgeStrechRel';...
    'SumIntAnglesFol';'SumIntAnglesRel';'SumExtAnglesFol';'SumExtAnglesRel' };
writeHeader(fileMassDist, headersMassDist);

fileAngles = strcat(folderEnergy, '/','Angles.csv');
if exist(fileAngles, 'file')
    delete(fileAngles) % always start with new file
end
headersAngles = {'HingeNumber'; 'All Angles'};
writeHeader(fileAngles, headersAngles);

fileMetadata = strcat(folderEnergy, '/','metadata.txt');
if exist(fileMetadata, 'file')
    delete(fileMetadata) % always start with new file
end
copyfile([folderResults '/metadata.txt'],fileMetadata);

function writeHeader(file, headers)

fid = fopen(file, 'wt') ;                         % Opens file for writing.
for j = 1:size(headers,1)
    fprintf(fid, '%s', headers{j});
    if j ~= size(headers,1)
        fprintf(fid, ',');
    else
        fprintf(fid, '\n');
    end
end
fclose(fid) ;                                          % Closes file.


function [hingeSet, opening] = getHingeSet(fileName)

parsedName = strsplit(fileName(1:end-4), '_');
hingeSetStr = parsedName{1};
if length(hingeSetStr)>2
    if strcmp(hingeSetStr(end-1:end),'op')
        hingeSetStr = hingeSetStr(1:end-2);
        opening = 1;
    else
        opening = 0;
    end
else
    opening = 0;
end
hingeSetStr = strrep(hingeSetStr, '[', '');
hingeSetStr = strrep(hingeSetStr, ']', '');
hingeSetStr = strsplit(hingeSetStr(1:end), ' ');
hingeSet = str2double(hingeSetStr)';


function [CM, Radios, Stdev, EhingeInt, SumIntAngles, SumExtAngles, maxStrech, minStrech] =...
    getData(extrudedUnitCell, opt, result)
CM = zeros(2,2,3);
Radios = zeros(2,2);
Stdev = zeros(2,2);
EhingeInt = zeros(2,2);
SumIntAngles = zeros(2,2);
SumExtAngles = zeros(2,2);
maxStrech = zeros(2,2);
minStrech = zeros(2,2);

foldingIterations = length(result.deform)-2;
for iter = 1:2
    currIter = foldingIterations+iter;
    startPos = extrudedUnitCell.node + result.deform(currIter).interV(1).V;
    endPos = extrudedUnitCell.node + result.deform(currIter).interV(end).V;
    CM(:,iter,:) = [mean(startPos); mean(endPos)]; %CM(inter, iter,:)
    startAllRad = sqrt(sum(abs(startPos-CM(1,iter)).^2,2));
    endAllRad = sqrt(sum(abs(endPos-CM(2,iter)).^2,2));
    Radios(:,iter) = [mean(startAllRad);mean(endAllRad)];
    Stdev(:,iter) = [std(startAllRad);std(endAllRad)];
    [EhingeInt(1,iter),SumIntAngles(1,iter), SumExtAngles(1,iter)] =...
        getIntEnergy(result.deform(currIter).interV(1), opt, extrudedUnitCell);
    [EhingeInt(2,iter),SumIntAngles(2,iter), SumExtAngles(2,iter)] =...
        getIntEnergy(result.deform(currIter).interV(end), opt, extrudedUnitCell);
    [maxStrech(1,iter), minStrech(1,iter)] = ...
        getExtremeStreching(result.deform(currIter).interV(1).Ve, opt, extrudedUnitCell);
    [maxStrech(2,iter), minStrech(2,iter)] = ...
        getExtremeStreching(result.deform(currIter).interV(end).Ve, opt, extrudedUnitCell);
%     end
end

function [EhingeIntSum, suminttheta, sumexttheta] = getIntEnergy(result, opt, extrudedUnitCell)
theta=result.theta;
EhingeInt = zeros(size(extrudedUnitCell.innerHinges,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[result.Ve(1:3:end) result.Ve(2:3:end) result.Ve(3:3:end)];

for i= 1:length(extrudedUnitCell.innerHinges)
    EhingeInt(i)=1/2*opt.Khinge*(theta(extrudedUnitCell.innerHinges(i))-...
        extrudedUnitCell.theta(extrudedUnitCell.innerHinges(i)))^2;
end

EhingeIntSum = sum(EhingeInt);
suminttheta = sum(theta(extrudedUnitCell.innerHinges));
theta(extrudedUnitCell.innerHinges) = [];
sumexttheta = sum(theta);

function [maxStrech, minStrech] = getExtremeStreching(u, opt, extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=(L-extrudedUnitCell.edgeL(i))/extrudedUnitCell.edgeL(i);            
end

maxStrech = max(dEdge);
minStrech = min(dEdge);

function normal = getavNormal(nodes,prevnormal) %%%%This can be used to get the normal of the faces

center = mean(nodes,1);
node1 = 1;
node2 = ceil(size(nodes,1)/2);

a=nodes(node1,:)-center;
b=nodes(node2,:)-center;
alpha=acos(sum(a.*b)/(norm(a)*norm(b)));
if imag(alpha) > 0
    alpha = 0;
end
normal=cross(a,b)/(norm(a)*norm(b)*sin(alpha));
if sum(isinf(normal))
    normal = prevnormal;
end

