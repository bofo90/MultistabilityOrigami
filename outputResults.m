function outputResults(unitCell,extrudedUnitCell,result,opt,filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT INFO FOR IMPLEMENTATION OF NEW GEOMETRIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col=hsv(32);
colt=col([12,26,31,20],:);
colt(1:2,:)=colt(1:2,:)+(1-colt(1:2,:))*0.3;
colt(3:4,:)=colt(3:4,:)+(1-colt(3:4,:))*0.1;
colt=[colt; col(2,:)+(1-col(2,:))*0.85];

colt(4,:)=[0,153,255]/255;
colt(5,:)=[225,225,225]/255;

%col=hsv(32);
%col=col+(1-col)*0;
%colt=col([32,31,8,13],:)
%colt=[colt; col(2,:)+(1-col(2,:))*0.85];
    
viewCoor=[sind(opt.AZ) -cosd(opt.AZ) sind(opt.EL)];
opt.tranPol=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARE PLOTTING UNDEFORMED CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if output folder is required, and create it if it doesn't exist
nameFolder=[pwd,'/Results/',opt.template,'/',opt.saveFile,'/images'];
if exist(nameFolder, 'dir')==0
    mkdir(nameFolder)
end

nref=size(unitCell.l,1);
if nref==0
    extrudedUnitCell.ref=[];
end

%prepare lattice vectors for plotting periodic structures and lattice
for ne=1:length(unitCell.Polyhedron)
    unitCell.Polyhedron(ne).latVec=detLatVec(unitCell.lhat,opt);
    unitCell.PolyhedronNew(ne)=unitCell.Polyhedron(ne);
end
extrudedUnitCell.latVec=detLatVec(unitCell.l,opt); 
%Prepare polyhdra and solid faces for efficient plotting
for ne=1:length(unitCell.Polyhedron)
    plotunitCell.Init(ne)=prepEffPlot(unitCell.Polyhedron(ne),viewCoor);
    plotunitCell.InitNew(ne)=prepEffPlot(unitCell.PolyhedronNew(ne),viewCoor);
end
plotextrudedUnitCell=prepEffPlot(extrudedUnitCell,viewCoor);

%Plot solid face with 100% transparency
f=figure('Position', [0 0 800 800]); hold on
for nc=1:size(extrudedUnitCell.latVec,1)
    for i=3:15
        c=(plotextrudedUnitCell.polFace(i).normal*viewCoor')>0;
        hs{nc,i}=patch('Faces',plotextrudedUnitCell.polFace(i).nod,'Vertices',plotextrudedUnitCell.lat(nc).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1.0,'edgealpha',1.0);
    end
end

%Plot polyhedra with 100% transparency
for ne=1:length(unitCell.Polyhedron)
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for i=3:15
            [c, plotunitCell.Init(ne).polFace(i).indexe, b1]=intersect(plotunitCell.Init(ne).polFace(i).index,unitCell.Polyhedron(ne).extrude);
            [c, plotunitCell.Init(ne).polFace(i).indexs, b1]=intersect(plotunitCell.Init(ne).polFace(i).index,unitCell.Polyhedron(ne).solidify);
            hie{ne,nc,i}=patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexe,:),'Vertices',plotunitCell.Init(ne).lat(nc).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',0.0,'edgealpha',0.0);
            his{ne,nc,i}=patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexs,:),'Vertices',plotunitCell.Init(ne).lat(nc).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',0.0,'edgealpha',0.0);
        end
    end
end
%Set axis
axis tight
xlim=1.1*get(gca,'xlim');
ylim=1.1*get(gca,'ylim');
zlim=1.1*get(gca,'zlim');
set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
hl2=plotOpt(opt);

opt.xlim=xlim;
opt.ylim=ylim;
opt.zlim=zlim;
if strcmp(opt.analysis,'result') || strcmp(opt.analysis,'savedata') || strcmp(opt.analysis,'plot')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PREPARE PLOTTING OF DEFORMED CONFIGURATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if result.numMode>0
        %prepare reference vectors for plotting periodic structures  
        refVec=detrefVec(unitCell,result,opt,extrudedUnitCell.ref);

        %Determine node coordinates and orientation of faces for each step in
        %the mode
        maxAxis=[0 0 0];
        minAxis=[0 0 0];
        Imax=2;
        Imin=-1;
        for nMode=1:result.numMode
            fram=0;
            for framMode=[0:size(result.deform(nMode).interV,2)-1]
                fram=fram+1;
                for nc=1:size(extrudedUnitCell.latVec,1)
                    plotextrudedUnitCell.mode(nMode).frame(fram).lat(nc).coor=extrudedUnitCell.node+result.deform(nMode).interV(framMode+1).V+ones(size(extrudedUnitCell.node,1),1)*(extrudedUnitCell.latVec(nc,:)-refVec(nMode).interV(framMode+1).val(nc,:));
                    ma=max(plotextrudedUnitCell.mode(nMode).frame(framMode+1).lat(nc).coor);
                    mi=min(plotextrudedUnitCell.mode(nMode).frame(framMode+1).lat(nc).coor);
                    maxAxis(maxAxis<ma)=ma(maxAxis<ma);
                    minAxis(minAxis>mi)=mi(minAxis>mi);
                end
                for i=3:15
                    plotextrudedUnitCell.mode(nMode).frame(fram).polFace(i).normal(1,:)=[0 0 0];
                    for j=1:size(plotextrudedUnitCell.polFace(i).nod,1)
                        n1a=plotextrudedUnitCell.mode(nMode).frame(fram).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,2),1:3)-plotextrudedUnitCell.mode(nMode).frame(framMode+1).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,1),1:3);
                        n2a=plotextrudedUnitCell.mode(nMode).frame(fram).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,3),1:3)-plotextrudedUnitCell.mode(nMode).frame(framMode+1).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,1),1:3);
                        n3a=cross(n1a,n2a);
                        plotextrudedUnitCell.mode(nMode).frame(fram).polFace(i).normal(j,:)=n3a;
                    end
                end
            end
        end
        %Update axis
        xlim=[minAxis(1),maxAxis(1)];
        ylim=[minAxis(2),maxAxis(2)];
        zlim=[minAxis(3),maxAxis(3)];
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
    end
end

if strcmp(opt.analysis,'info')
    %First make solid face with 100% transparency
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:15
            set(hs{nc,i},'facealpha',0,'edgealpha',0);
        end
    end
    %INTERNAL POLYHEDRA TESSELATION
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for ne=1:length(unitCell.Polyhedron)
            for i=3:15       
                set(hie{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                set(his{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
            end               
        end
    end
    printHigRes(f,opt,[opt.template,'_internalPolyhedra'],nameFolder); 
    %EXTRUDED TESSELATION
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:15
            set(hs{nc,i},'facealpha',1,'edgealpha',1)         
        end
    end
    for ne=1:length(unitCell.Polyhedron)
        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            for i=3:15
                set(hie{ne,nc,i},'facealpha',0,'edgealpha',0);
                set(his{ne,nc,i},'facealpha',0,'edgealpha',0);
            end
        end
    end
    printHigRes(f,opt,[opt.template,'_extrudedMat'],nameFolder);
    %PLOT ALL INFORMATION UNIT CELL
    f=figure('Position', [0 0 800 800]);
    hold on;
    for ne=1:length(unitCell.Polyhedron)
        for i=3:15
            patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexe,:),'Vertices',plotunitCell.Init(ne).lat(1).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
            patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexs,:),'Vertices',plotunitCell.Init(ne).lat(1).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
        end
    end
     
    for ne=1:length(unitCell.Polyhedron)
        coorCenter=sum(unitCell.Polyhedron(ne).node)/size(unitCell.Polyhedron(ne).node,1);
        text(coorCenter(1),coorCenter(2),coorCenter(3),num2str(ne),'fontsize',30,'color','g')
        for i=1:size(unitCell.Polyhedron(ne).node,1)
            coor=unitCell.Polyhedron(ne).node(i,:);
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*k')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20)
        end
        for i=1:size(unitCell.Polyhedron(ne).face,1)
            coor=sum(unitCell.Polyhedron(ne).node(unitCell.Polyhedron(ne).face{i},:))/length(unitCell.Polyhedron(ne).face{i});
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*b')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20,'color','b')
        end
    end 
    hl2=plotOpt(opt); 
    axis tight
    printHigRes(f,opt,[opt.template,'_internalPolyhedraVertFaces'],nameFolder); 
    %EXTRUDED STATE UNIT CELL
    f4=figure('Position', [0 0 800 800]);
    hold on;
    for i=3:15
        c=(plotextrudedUnitCell.polFace(i).normal*viewCoor')>0;
        hsu{1,i} = patch('Faces',plotextrudedUnitCell.polFace(i).nod,'Vertices',plotextrudedUnitCell.lat(1).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1,'edgealpha',1);
    end

    hl2=plotOpt(opt); 
    axis tight
    printHigRes(f4,opt,[opt.template,'_extrudedPolyhedra'],nameFolder); 
    %EXTRUDED STATE 'INFO' VERTEX NUMBERS
    for i=3:15
        set(hsu{1,i},'facealpha',0.75,'edgealpha',0.75)         
    end
    coorCenter=sum(extrudedUnitCell.node)/size(extrudedUnitCell.node,1);    
    for i=1:size(extrudedUnitCell.node,1)
            coor=extrudedUnitCell.node(i,:);
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            hline{i} = line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':');
            hmark{i} = plot3(coor(1),coor(2),coor(3),'*k');
            htext{i} = text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20);
    end
    printHigRes(f4,opt,[opt.template,'_extrudedPolyhedraVert'],nameFolder); 
    %EXTRUDED STATE 'INFO' ANGLE NUMBERS
    for i=1:size(extrudedUnitCell.node,1)
        delete(hline{i});
        delete(hmark{i});
        delete(htext{i});
    end
    for i=1:size(extrudedUnitCell.nodeHingeEx,1)
            coor1=extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,1),:);
            coor2=extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,2),:);
            coor=coor1/2+coor2/2;
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*k')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20)
    end
    printHigRes(f4,opt,[opt.template,'_extrudedPolyhedraEdges'],nameFolder);


end

if strcmp(opt.analysis,'plot')
        %First make solid face with 100% transparency
        for nc=1:size(extrudedUnitCell.latVec,1)
            for i=3:15
                set(hs{nc,i},'facealpha',0,'edgealpha',0);
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT SPACE-FILLING ASSEMBLY POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            fram=fram+1;
            for ne=1:length(unitCell.Polyhedron)
                for i=3:15
                    set(hie{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                    set(his{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                end
                set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim); 
            end
        end
        axis equal
%         printHigRes(f,opt,'polyhedra',nameFolder)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT MOVING OF POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %New coordinates

        for ne=1:length(unitCell.Polyhedron)
            for nc=1:size(unitCell.Polyhedron(ne).latVec,1)
                plotunitCell.Init(ne).lat(nc).coorNew(:,1)=plotunitCell.InitNew(ne).lat(nc).coor(:,1)+extrudedUnitCell.latVec(nc,1)-unitCell.PolyhedronNew(ne).latVec(nc,1);
                plotunitCell.Init(ne).lat(nc).coorNew(:,2)=plotunitCell.InitNew(ne).lat(nc).coor(:,2)+extrudedUnitCell.latVec(nc,2)-unitCell.PolyhedronNew(ne).latVec(nc,2);
                plotunitCell.Init(ne).lat(nc).coorNew(:,3)=plotunitCell.InitNew(ne).lat(nc).coor(:,3)+extrudedUnitCell.latVec(nc,3)-unitCell.PolyhedronNew(ne).latVec(nc,3);
                        %Update position
                for i=3:15
%                             set(hie{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
%                             set(his{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
                    set(hie{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
                    set(his{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
                end 
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
%         printHigRes(f,opt,'Polyhedra_Packing_Expanded',nameFolder)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT SELECTED FACES TO EXTRUDE, SOLIDIFY AND REMOVE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            for ne=1:length(unitCell.Polyhedron)
                for i=3:15
                    if ~isempty(plotunitCell.Init(1).polFace(i).indexs)
                        c=(plotunitCell.Init(ne).polFace(i).normal(plotunitCell.Init(ne).polFace(i).indexs,:)*viewCoor')>0;
                        set(his{ne,nc,i},'facealpha',1,'edgealpha',1,'facevertexCData',(c*colt(4,:)+abs(1-c)*colt(5,:)));
                    end
                end
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);  
%         printHigRes(f,opt,'facetype',nameFolder)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT EXTRUSION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nc=1:size(extrudedUnitCell.latVec,1)
            for i=3:15
                set(hs{nc,i},'facealpha',1,'edgealpha',1)      
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %REMOVE POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ne=1:length(unitCell.Polyhedron)
            for nc=1:size(unitCell.Polyhedron(1).latVec,1)
                for i=3:15
                    set(hie{ne,nc,i},'facealpha',0,'edgealpha',0);
                    set(his{ne,nc,i},'facealpha',0,'edgealpha',0);
                end
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
        printHigRes(f,opt,[filename,'_0_undeformed'],nameFolder);     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT MODES INDIVIDUALLY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nMode=1:result.numMode
            pause(1)
%             fprintf('end part');
            if strcmp(opt.gethistory,'on')
                for framMode=1:length(plotextrudedUnitCell.mode(nMode).frame)
    %                 framMode = length(plotextrudedUnitCell.mode(nMode).frame);
                    for nc=1:size(extrudedUnitCell.latVec,1)
                        for i=3:15    
                            c=(plotextrudedUnitCell.mode(nMode).frame(framMode).polFace(i).normal*viewCoor'>0);
                            viewCoor=get(gca,'view');
                            viewCoor=[sind(viewCoor(1)) -cosd(viewCoor(1)) sind(viewCoor(2))];
                            set(hs{nc,i},'Vertices',plotextrudedUnitCell.mode(nMode).frame(framMode).lat(nc).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1.0);
                        end
                    end
                    printGif(opt,framMode,f,nameFolder,[filename,'_',num2str(nMode),'_deformed']);%'_',sprintf('%2.3f_%2.3f_%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge),
                    if framMode==length(plotextrudedUnitCell.mode(nMode).frame)
                        printHigRes(f,opt,[filename,'_',num2str(nMode),'_deformed'],nameFolder);%'_',sprintf('%2.3f_%2.3f_%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge),
                    end
                end
            else
                framMode = length(plotextrudedUnitCell.mode(nMode).frame);
                for nc=1:size(extrudedUnitCell.latVec,1)
                    for i=3:15     
                        c=(plotextrudedUnitCell.mode(nMode).frame(framMode).polFace(i).normal*viewCoor'>0);
                        viewCoor=get(gca,'view');
                        viewCoor=[sind(viewCoor(1)) -cosd(viewCoor(1)) sind(viewCoor(2))];
                        set(hs{nc,i},'Vertices',plotextrudedUnitCell.mode(nMode).frame(framMode).lat(nc).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1.0);
                    end
                end
                if framMode==length(plotextrudedUnitCell.mode(nMode).frame)
                    printHigRes(f,opt,[filename,'_',num2str(nMode),'_deformed'],nameFolder);%'_',sprintf('%2.3f_%2.3f_%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge),
                end
            end
            
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USED FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function latVec=detLatVec(perVec,opt)
    latVec=[];
    nref=size(perVec,1);
    switch nref
        case 0
            latVec=[0 0 0];                 
        case 1
            for i=1:opt.plotPer
                latVec=[latVec; (i-1)*perVec(1,:)];
            end
        case 2
            for j=1:opt.plotPer
                for i=1:opt.plotPer
                    latVec=[latVec; (i-1)*perVec(1,:)+(j-1)*perVec(2,:)];
                end
            end
        case 3
            for k=1:opt.plotPer
                for j=1:opt.plotPer
                    for i=1:opt.plotPer
                        latVec=[latVec; (i-1)*perVec(1,:)+(j-1)*perVec(2,:)+(k-1)*perVec(3,:)];
                    end
                end
            end     
    end  
    
function refVec=detrefVec(unitCell,result,opt,ref)
    for nMode=1:length(result.deform)
        for fram=0:size(result.deform(nMode).interV,2)-1
            refVec(nMode).interV(fram+1).val=[];
            switch size(unitCell.l,1);
                case 0
                    refVec(nMode).interV(fram+1).val=[0 0 0];                    
                case 1
                    for i=1:opt.plotPer
                        refVec(nMode).interV(fram+1).val=[refVec(nMode).interV(fram+1).val; (i-1)*result.deform(nMode).interV(fram+1).V(ref(1),:)];
                    end
                case 2
                    for j=1:opt.plotPer
                        for i=1:opt.plotPer
                            refVec(nMode).interV(fram+1).val=[refVec(nMode).interV(fram+1).val; (i-1)*result.deform(nMode).interV(fram+1).V(ref(1),:)+(j-1)*result.deform(nMode).interV(fram+1).V(ref(2),:)];
                        end
                    end
                case 3
                    for k=1:opt.plotPer
                        for j=1:opt.plotPer
                            for i=1:opt.plotPer
                                refVec(nMode).interV(fram+1).val=[refVec(nMode).interV(fram+1).val; (i-1)*result.deform(nMode).interV(fram+1).V(ref(1),:)+(j-1)*result.deform(nMode).interV(fram+1).V(ref(2),:)+(k-1)*result.deform(nMode).interV(fram+1).V(ref(3),:)];
                            end
                        end
                    end     
            end 
        end
            
    end
    
function plotg=prepEffPlot(som,viewCoor)
    for i=3:20
        plotg.polFace(i).nod=[];
        plotg.polFace(i).index=[];
    end
    for i=1:length(som.face)
        l=length(som.face{i});
        plotg.polFace(l).nod=[plotg.polFace(l).nod;som.face{i}];
        plotg.polFace(l).index=[plotg.polFace(l).index; i];
    end

    %COORDINATES AND FACE NORMAL POLYHEDRON
    for nc=1:size(som.latVec,1)
        plotg.lat(nc).coor=som.node+ones(size(som.node,1),1)*(som.latVec(nc,:));
        if isfield(som,'nodeNew')
            plotg.lat(nc).coorNew=som.nodeNew+ones(size(som.nodeNew,1),1)*(som.latVec(nc,:));
        end
    end
    for i=3:20
        plotg.polFace(i).normal(1,:)=[0,0,0];
        for j=1:size(plotg.polFace(i).nod,1)
            n1a=plotg.lat(nc).coor(plotg.polFace(i).nod(j,2),1:3)-plotg.lat(nc).coor(plotg.polFace(i).nod(j,1),1:3);
            n2a=plotg.lat(nc).coor(plotg.polFace(i).nod(j,3),1:3)-plotg.lat(nc).coor(plotg.polFace(i).nod(j,1),1:3);
            n3a=cross(n1a,n2a);
            plotg.polFace(i).normal(j,:)=n3a;
        end
    end    
    
function hl2=plotOpt(opt)
    lighting flat; 
    axis equal; 
    axis off; 
    view(opt.AZ,opt.EL);
    hl2=camlight(-15,30); 
    set(gcf,'color','w');
    hl2.Style='infinite';
    material dull;
    set(gca,'cameraviewanglemode','manual');
    %set(gca,'outerposition',[0 0 1.0 1.0])
    
function printGif(opt,fram,f,nameFolder,nam)
    pause(1/opt.frames)
    name=[nameFolder,'/', nam];
    switch opt.analysis
        case 'modes'
            name=[name,'_Modes_'];
    end
    switch opt.periodic
        case 'on'
            name=[name,'_pcb_',num2str(opt.plotPer),'.gif'];
        case 'off'
            name=[name,'.gif'];
    end
    frame = getframe(f.Number);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if fram==1
        imwrite(imind,cm,name,'gif', 'Loopcount',inf,'delaytime',0.04);
    else
        imwrite(imind,cm,name,'gif','WriteMode','append','delaytime',0.04);
    end
    
function printHigRes(f,opt,nam,nameFolder)
    pause(1/opt.frames)
    name=[nameFolder,'/',nam,'.png'];
    savefig([nameFolder,'/',nam])
    figpos=getpixelposition(f); %dont need to change anything here
    resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,...
    'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
    print(f,name,'-dpng',['-r',num2str(opt.figDPI)],'-opengl') %save file

function [f,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his)
    f=figure('Position', [0 0 800 800]);
    ax = axes;
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:15
            hs{nc,i}=copyobj(hs{nc,i},ax);
        end
    end
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for ne=1:length(unitCell.Polyhedron)
            for i=3:15
                hie{ne,nc,i}=copyobj(hie{ne,nc,i},ax);
                his{ne,nc,i}=copyobj(his{ne,nc,i},ax);
            end
        end
    end
%     set(gca,'xlim',opt.xlim,'ylim',opt.ylim,'zlim',opt.zlim);
    hl2=plotOpt(opt);