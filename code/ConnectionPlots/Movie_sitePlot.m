function WN_sitePlot(subj)
   if nargin==0
       subj='MG51';
   end
%    close all
   ih=1;
   h={'lh'};
   addpath('/usr/pubsw/packages/freesurfer/RH4-x86_64-R600/matlab')%read_curv for sulc
   addpath('/usr/pubsw/packages/MMPS/MMPS_235/matlab/fsurf_tools')
   addpath('/home/pubsw/packages/MMPS/MMPS_235/matlab/mmil_utils')
    [srfFSavg.(h{ih}).vertices, srfFSavg.(h{ih}).faces] = read_surf(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/electrode_location_files/fsaverage/surf/%s.pial',h{ih}));
    srfFSavg.(h{ih}).faces=srfFSavg.(h{ih}).faces+1;
%     srfFSavg.(h{ih}).vertNorm = patchnormals(srfFSavg.(h{ih}));
    annotfil2 = sprintf('~/MULTI/fsaverage/%s.HCP-MMP1.annot',h{ih});
    % [~,Lparc,parc]=read_annotation(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/electrode_location_files/fsaverage/label/%s.aparc.a2009s.annot',h{ih}));
    [~,LparcFSavg,parcFSavg]=read_annotation(annotfil2);
    [~,LparcFSavg]=ismember(LparcFSavg,parcFSavg.table(:,5));
    
    nParc=parcFSavg.numEntries;
    
    tmpTab=[1 1 1; parcFSavg.table(:,1:3)];
    val{ih}=tmpTab(LparcFSavg+1,:);
%     val{ih}(all(val{ih}==0,2),:) = 128;
%     val{ih} = val{ih} - (val{ih} - mean(val{ih},2))/(4/3);
    
    %mobj = matfile('/home/bqrosen/projects/SS_SIM//out/fsaverage_ico7_pial_surf_lh.mat');
    %mobj = matfile(sprintf('/home/bqrosen/projects/SS_SIM//out/fsaverage_ico7_inflated_surf_%s.mat',h{ih}));

    h1=figure('Position',get(0,'Screensize'));
%     srfFSavg.(h{ih}).faces=srfFSavg.(h{ih}).faces+1;
%     srf.(h{ih}).vertNorm = patchnormals(srf.(h{ih}));
    %pH = triplot(srf.(h{ih}).vertices,srf.(h{ih}).faces,val{ih}(:,1))%,'faces_skin');
    pH = trisurf(srfFSavg.(h{ih}).faces,srfFSavg.(h{ih}).vertices(:,1),srfFSavg.(h{ih}).vertices(:,2),srfFSavg.(h{ih}).vertices(:,3),...
        'LineStyle','none','FaceVertexCData',val{ih}/255,'CDataMapping','direct');
    %xlim([min(pH.XData(:)) max(pH.XData(:))]);ylim([min(pH.YData(:)) max(pH.YData(:))]);zlim([min(pH.ZData(:)) max(pH.ZData(:))])
    axis image off
    % figure background color
    set(gcf,'color',[0.5 0.5 0.5])
    %myaa(9)
    
    % patch properties
    %pH.FaceColor = [.5 .5 .5];
    pH.FaceAlpha = 1;
    hold on
    view(90+180*(ih-1),0)
    delete(findobj(gcf,'type','light'));
    camlight headlight;material dull;lighting gouraud
    view(180*(ih-1)-90,0);camlight headlight
    title([subj ' ' h{ih}])
    
%     distFun=@(a,b) sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2 + (a(3)-b(3))^2);
    
    load(sprintf('/home/jgarret/TaskAnalysis/lang/MG/parcs/%s_%s_iEEG_HCP_fsAvg.mat',subj,h{ih}),'eegParcHCP','chanCoordFSavg')
    vertCtr=mean(srfFSavg.(h{ih}).vertices,1);
    
    for n = 1:size(eegParcHCP,1)

        %figure(figNum+sumTab.hemi(n))

%         eleVert=find(annot.Lidx.(h{ih})==n,1,'first');
%         vertCoords=srf.(h{ih}).vertices(eleVert,:) - srf.(h{ih}).vertNorm(eleVert,:)/10;
        vertCoords=cell2mat(chanCoordFSavg(n,2:4));

    %     vertCoords = srf.(h{sumTab.hemi(n)}).vertices(sumTab.icoVert(n),:)...
    %         - srf.(h{sumTab.hemi(n)}).vertNorm(sumTab.icoVert(n),:);

%         if strcmp(subj,'NY240v2')
%             vertCoords=siteCoords(n,:);% - srf.(h{ih}).vertNorm(nearSubjVert(c),:)/10;
            
%         end
        
        clr=[0.75 0 1];
        
            vertCoords = vertCoords + 3*(vertCoords-vertCtr)/norm(vertCoords-vertCtr);% + srf.(h{ih}).vertNorm(nearSubjVert(c),:)/20;
        plot3(vertCoords(1),vertCoords(2),vertCoords(3),'.','markersize',45,'color',[1 0.5 0]);
        if ~isempty(vertCoords)
            text(vertCoords(1),vertCoords(2),vertCoords(3),eegParcHCP{n,2},'FontSize',10,'Color',clr,'FontWeight','bold');
           
        end


    end
    
    
end
