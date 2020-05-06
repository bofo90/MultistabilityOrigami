function selectHinges(unitCell, extrudedUnitCell, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  Carefull with the selection of hinges! There is a problem when the
%%%%  amount of hinges per vertex is the same as the amount of hinges per
%%%%  face if all faces are equal in the same vertex (1 flavour in ). For
%%%%  example the tetrahedron (3 hinges per vertix and each face has 3
%%%%  hinges). There are specific hinge selections that from the distance
%%%%  matriz have the same eigenvalu but there are different.
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the function for selecting hinges
if strcmp(opt.analysis, 'selecthinges')
    fprintf('Creating graph...\n');
    [G] = buildGraph(unitCell, extrudedUnitCell);
    fprintf('Calculating distances...\n');
    dis = getDistance(G);
    % first generate the list of all hinge sets
    fprintf('Writting hinges...\n');
    getAllHinges(G, dis, opt);
    fprintf('Ploting hinges...\n');
    if strcmp(opt.createFig, 'on')
        plotHinges(unitCell, extrudedUnitCell, opt);
    end
    fprintf('Done with selecting hinges...\n');
end

% hi all
function [G] = buildGraph(unitCell, extrudedUnitCell)

numEdges = size(unitCell.Polyhedron.edge, 1);
G_adjacency = zeros(numEdges);

hingeNames = double.empty(numEdges, 0); % corresponds to names in extrudedUnitCell

flavourlist = double.empty(0,2);
hingeTypes = double.empty(numEdges, 0);
orderedEdges = sort(unitCell.Polyhedron.edge,2);
orderedFaces = sortFace(unitCell.Polyhedron);

for ii = 1:numEdges
    n1 = unitCell.Polyhedron.edge(ii,1);
    n2 = unitCell.Polyhedron.edge(ii,2); % endnodes of the edge
    % getting the edge numbering from nodeHingeEx (consistent with how
    % they're numbered in outputResults
    [~,idx1]=ismember(orderedEdges(ii,:),extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
%     [~,idx2]=ismember([n2 n1],extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
    hingeNames(ii) = idx1;%+idx2;
    
    % getting the types of faces associated with current edge
    faceTypes = [0 0];
    nextEdge = [0 0;0 0];
    for ff = 1:length(unitCell.Polyhedron.face)
        n_face = orderedFaces{ff};
        pos1 = find(n1==n_face);
        pos2 = find(n2==n_face);
        sides = length(n_face);
        if sum(pos1) && sum(pos2)
            index = find(0==faceTypes, 1);
            faceTypes(index) = sides;
            
            maxPos = max(pos1,pos2);
            if maxPos == sides
                if min(pos1,pos2) == sides-1
                    nextEdge(index,:) = [n_face(end) n_face(1)];
                else
                    nextEdge(index,:) = [n_face(1) n_face(2)];
                end
            else
                nextEdge(index,:) = [n_face(maxPos) n_face(maxPos+1)];
            end

        end
    end
    [added, flavour] = ismember([min(faceTypes) max(faceTypes)],flavourlist, 'rows');
    if ~added
       flavourlist(end+1,:) = [min(faceTypes) max(faceTypes)];
       flavour =  size(flavourlist,1);
    end
    hingeTypes(ii) = [flavour];
    
    nextEdge = sort(nextEdge,2);
    [~,adedge1] = ismember(nextEdge(1,:),orderedEdges,'rows');
    [~,adedge2] = ismember(nextEdge(2,:),orderedEdges,'rows');
    G_adjacency(ii, adedge1) = 1;
    G_adjacency(ii, adedge2) = 1;
end

G_adjacency = sparse(G_adjacency);

% finally, make a graph of the hinges
G = digraph(G_adjacency);

% assign names and other properties to the nodes of the graph (edges on the polyherdon)
G.Nodes.Names = hingeNames';
G.Nodes.Type = hingeTypes';

function orderedFaces = sortFace(Polyhedron)

center = mean(Polyhedron.node, 1);
orderedFaces = Polyhedron.face;

for i = 1:size(Polyhedron.face,1)
    a = Polyhedron.node(Polyhedron.face{i}(1),:)-Polyhedron.node(Polyhedron.face{i}(2),:);
    b = Polyhedron.node(Polyhedron.face{i}(3),:)-Polyhedron.node(Polyhedron.face{i}(2),:);
    node = Polyhedron.node(Polyhedron.face{i}(2),:)-center;
    direction = dot(node, cross(a,b));
    if direction < 0
        orderedFaces{i} = flip(orderedFaces{i});
    end
end


function dis = getDistance(G)
% Returns the ``distance'' matrix with flavours
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - an digraph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% dis - a distance matrix
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% initialising...
H = height(G.Nodes);
dis = zeros(H);

% loop through everything, because we need to
for ii = 1:H
    for jj = 1:H
        % get all possible paths between two nodes *ii* and *jj* both ways
        pathsiijj = pathbetweennodes(G, ii, jj);
      
        % get the ``distance'' from node *ii* to node *jj*
        % by looping through all minimum length paths
        pathDis = inf; % initialise with infinity
        for mCt = 1:size(pathsiijj,1)
            currentPath = pathsiijj(mCt,:);
            currentPathTypes = G.Nodes.Type(currentPath);
            currentPathDisStr = strjoin(string(currentPathTypes),'');
            
            % convert to number and take the smallest value of the
            % ``distance'' and the flipped ``distance''
            currentPathDis = min(str2double(currentPathDisStr), str2double(reverse(currentPathDisStr)));
            if currentPathDis < pathDis
                pathDis = currentPathDis;
            end
        end

        % update distance matrix
        dis(ii,jj) = pathDis;
%         dis(jj,ii) = pathDis;
    end
end


function pth = pathbetweennodes(G, src, snk)
%PATHBETWEENNODES Return all paths between two nodes of a graph
%
% pth = pathbetweennodes(adj, src, snk)
% pth = pathbetweennodes(adj, src, snk, vflag)
%
%
% This function returns all simple paths (i.e. no cycles) between two nodes
% in a graph.  Not sure this is the most efficient algorithm, but it seems
% to work quickly for small graphs, and isn't too terrible for graphs with
% ~50 nodes.
%
% Input variables:
%
%   adj:    adjacency matrix
%
%   src:    index of starting node
%
%   snk:    index of target node
%
%   vflag:  logical scalar for verbose mode.  If true, prints paths to
%           screen as it traverses them (can be useful for larger,
%           time-consuming graphs). [false]
%
% Output variables:
%
%   pth:    cell array, with each cell holding the indices of a unique path
%           of nodes from src to snk.

% input validation added by yun, May 02, 2017
if src==snk
    pth = [src];
    return;
end
pth = 1:size(G.Nodes,1);
for node = pth
    next(node,:) = successors(G,node)';
end 
%Recursive search of shortest distance between src and snk
[pth,~,~,~] = nextNode(pth, src, snk, next);
%Recursive search of shortest distance between snk and src
% [pth,~,~,~] = nextNode(pth, snk, src, next);


function [paths, stack, snk, next] = nextNode(paths, stack, snk, next)

%Check if the distance is longer than the previous found one
if length(stack) > length(paths(end,:))
    return;
end
%Check if we dont have loops in the stack
if any(diff(sort(stack))==0)
    return;
end
%Check if you already end in the finnal edge
if stack(end) == snk
    %if the distance is shorter, erease the rest, if not just stack them
    if length(stack) < length(paths(end,:)) 
        paths = stack;
    else
        paths = [paths; stack];
    end
    return;
end
%go through the next edges and call the function again
originalstack = stack;
for j = next(originalstack(end),:)
    newstack = [originalstack j];
    [paths, stack, snk, next] = nextNode(paths, newstack, snk, next);
end

function getAllHinges(G, dis, opt)
% get the complete list of hinge sets, and dump them in a csv file
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - a directed graph created from the options
% flavourTypes
% flavourNum
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% A .csv file of the complete list of hinge sets, where each row is one
% such hinge set
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun

% name of all hinges converted to an array of numbers
allHinges = 1:height(G.Nodes);

% starts generating all hignes
hingeSetsPrev(height(G.Nodes)+1).all = allHinges;
for N = 1:height(G.Nodes)-1  %height(G.Nodes)-1
    fprintf('Calculating distance for %d nodes\n', N);
    if N <= 0.5 * height(G.Nodes)
        hingeSets = chooseHinges(G, dis, N, hingeSetsPrev(N).all);
        % getting ready for the next loop
    end
%     if N == 0.5 * height(G.Nodes)
%         hingeSets = chooseHinges(G, dis, N, hingeSetsPrev(N).all);
%         for ct = 1:size(hingeSets, 1)
%             % write complement selections to file
%             conjSet = allHinges;
%             conjSet(hingeSets(ct,:)) = [];
%             hingeSets(end+1,:) = conjSet;
%         end
%     end
    if N > 0.5 * height(G.Nodes)
        Nconj = height(G.Nodes)-N;
        hingeSets = zeros(size(hingeSetsPrev(Nconj+1).all, 1), N);
        for ct = 1:size(hingeSets, 1)
            % write complement selections to file
            conjSet = allHinges;
            conjSet(hingeSetsPrev(Nconj+1).all(ct,:)) = [];
            hingeSets(ct,:) = conjSet;
        end
    end
    
    hingeSetsPrev(N+1).all = hingeSets;
    
    if N >= opt.maxHinges
        break;
    end
end

saveHinges(hingeSetsPrev, opt, G)

function saveHinges(hingeSetsPrev, opt, G)

% create folder if it doesn't exist
folderName = strcat(pwd, '/Results/hingeList_reduced/');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, opt.template, '.csv');
if exist(fileName, 'file')
    delete(fileName) % always start with new file
end

for i = 2:height(G.Nodes)+1
    if size(hingeSetsPrev(i).all,1) > 1
        dlmwrite(fileName, G.Nodes.Names(hingeSetsPrev(i).all), 'delimiter', ',', '-append');
    else
        dlmwrite(fileName, G.Nodes.Names(hingeSetsPrev(i).all)', 'delimiter', ',', '-append');
    end
end


function [hinges] = chooseHinges(G, dis, N, hingeSetsPrev)
% select *N* number of hinges from graph *G* to actuate, the flavoured
% version
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - graph object
% dis - matrix containing the minimum distance between any two nodes
% N - number of hinges to be actuated
% hingeSetsPrev - the complete hinge sets where each set contains (N-1) hinges
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hinges - an x*N array, where each row is a set of hinges to be actuated
%          at the same time
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun

% input validation
if N < 0 || N > size(G.Nodes,1)/2
    hinges = [];
    warning('You are selecting the wrong number of hinges to actuate!')
elseif N == 0
    hinges = [];
elseif N == 1
    % actuate different types of hinges
    [~,hinges,~] = unique(G.Nodes.Type);
else% N >= 2
    % initialising...
    numberHinges = size(G.Nodes,1);
    hinges = zeros(size(hingeSetsPrev, 1)*(numberHinges-N),N);
    distanceEig = zeros(size(hingeSetsPrev, 1)*(numberHinges-N),N);
    index = 1;
    % loop through each hinge set
    for ct = 1:size(hingeSetsPrev, 1)
        hingeSet = hingeSetsPrev(ct, :);
        candidateHinges = (hingeSet(end)+1):numberHinges;
        % then start comparing matrices
        for nodeCt = candidateHinges
            newSet = [hingeSet nodeCt];
            newDisMat = dis(newSet, newSet);
            newDisEig = round(1000*sort(eig(newDisMat))') / 1000; % round off 

            % if newDisEig is not in cache, add current solution
%             sameEig = sum(distanceEig(1:index,:) == newDisEig, 2);
            cachedistanceEig = distanceEig(1:index,:);
            for eigenvalue = 1:N
                sameEig = cachedistanceEig(:,eigenvalue) == newDisEig(eigenvalue);
                cachedistanceEig = cachedistanceEig(any(sameEig,2),:);
            end
            if size(cachedistanceEig,1) < 1

%             if ~sum(sameEig == N)
                hinges(index, :) = newSet;
                distanceEig(index, :) = newDisEig;
                index = index+1;
            end
        end
    end 
    
    hinges( ~any(hinges,2),:) = [];
end

function plotHinges(unitCell, extrudedUnitCell, opt)

folderName = strcat(pwd, '/Results/hingeList_reduced/');
if ~exist(folderName, 'dir')
    fprintf('There is no folder with hinge list');
    return;
end

fileName = strcat(folderName, opt.template, '.csv');
if ~exist(fileName, 'file')
    fprintf('There is no file for the current polyhedra with hinge list');
    return;
end

hingeList = csvread(fileName);

for ii = 1:size(hingeList,1)
    hinges = hingeList(ii, :);
    hinges = hinges(0~=hinges);
    
    f = plotSelectedHinges(hinges, unitCell, extrudedUnitCell, false);
    if strcmp(opt.saveFig, 'on')
        folderName = strcat(pwd, '/Results/hingeList_reduced/', opt.template, '/', ...
                     num2str(length(hinges)), '/');
        if ~exist(folderName, 'dir')
            mkdir(folderName);
        end
        saveas(f, strcat(folderName, mat2str(hinges), '.png'));
        
    end
    close(f)

end

function f = plotSelectedHinges(hingeList, unitCell, extrudedUnitCell, orderedNames)
% function plotSelectedHinges(geometry, hingeList)
% highlights the set of selected hinges in a given polyhedron, for better
% visualisation
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingelist - a list containing the index of the selected hinges
% unitCell
% extrudedUnitCell
% orderedNames - boolean indicating whether the hinge names are printed as
%                the original order as in 
%                extrudedUnitCell.nodeHingeEx --> orderedNames = false
%                or in ordered version --> orderedNames = true
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% f - a graph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on 30 Oct, 2017
% agustin


if ~exist('orderedNames', 'var')
    orderedNames = true;
end

vertices = unitCell.Polyhedron.node;

% plot transparent polyhedra with patch()
faces = unitCell.Polyhedron.face;
figure; hold on;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); view([20 10]); axis('equal')
axis('off')
xlabel('x'); ylabel('y'); zlabel('z')
if 1==length(hingeList)
    hingeStr = num2str(hingeList);
    hingeStr = strcat('[',hingeStr,']');
else
    hingeStr = mat2str(hingeList);
end

title(strcat('Hinges selected: ',hingeStr(2:(end-1))))
for ct = 1:length(faces)
    patch('Vertices',vertices, 'Faces',faces{ct},...
          'FaceColor',[.73 .86 .95], 'FaceAlpha',0.3, ...
          'EdgeColor',[.2 .2 .2], 'lineWidth',0.5);
end
light('Position',[1 -1 1],'Style','local')


% find corresponding nodes on the hingeList and highlight those edges
for ed = 1:length(hingeList)
    endNodes = extrudedUnitCell.nodeHingeEx(hingeList(ed), 1:2);
    X = vertices(endNodes, 1);
    Y = vertices(endNodes, 2);
    Z = vertices(endNodes, 3);
    line(X, Y, Z, 'lineStyle',':', 'color', [.01 .23 .45], 'LineWidth', 8)
end


% add lables for the edges
for ct = 1:size(unitCell.Polyhedron.edge, 1)
    endNodes = unitCell.Polyhedron.edge(ct, :);
    
    if orderedNames
        hingeName = num2str(ct);
    else
        [~, hinge] = ismember(endNodes,...
                     extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
        hingeName = num2str(hinge);
    end
    
    textPos = 0.5 * (vertices(endNodes(1), :) + vertices(endNodes(2), :));
    text(textPos(1)+.01, textPos(2)+.01, textPos(3)-.04, hingeName, ...
        'FontSize', 26, 'FontWeight', 'bold', 'color', [.2 .2 .2]);
end

f = gcf;
