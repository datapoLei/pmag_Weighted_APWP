clc, close all, clear all

% load data
[num,txt,raw]=xlsread('./Paleopoles_Pangea_Formation.xlsx',3);

% set projection: eqaconic, eqdcylin, ortho, mollweid, robinson, eqdazim
projtype='ortho'; viewA=[-15 15 0];

% set time steps (Tinv), running window size (window) and Smooth_Parameter
Tinv=10; window=30; Smooth_Parameter=1000;  

% flags for plotting results
meanEflag=0; % without APWP errors
meanEflag=1; % with APWP errors
AgeResNoflag=0; % index
AgeResNoflag=1; % Age
% AgeResNoflag=2; % ResNo
% AgeResNoflag=3; % No labels

%% Do not change the code below

data=cell2mat(raw(2:end,6:13)); data(1,:)=[];
ResNo=cell2mat(raw(2:end,1)); ResNo(1,:)=[];
poles=[data(:,4) data(:,1:3) ResNo data(:,8)]; poles=sortrows(poles,1);
Euler=[data(:,5:7) data(:,4)]; Euler=sortrows(Euler,4);

AgeAll=cell2mat(raw(2:end,3:4)); AgeAll(1,:)=[];
AgeAll=[AgeAll data(:,4)]; AgeAll=sortrows(AgeAll,3);

plateIDq=cell2mat(raw(2:end,5)); plateIDq(1,:)=[];
plateIDq=[plateIDq data(:,4)]; plateIDq=sortrows(plateIDq,2);

for i=1:length(poles(:,1))
    [lonMR,latMR]=step2_Wing_Rot(poles(i,3),poles(i,2),Euler(i,2),Euler(i,1),Euler(i,3));
    polesR(i,:)=[lonMR,latMR];
    
    % flip to south poles
    [latA, lonA]=antipode(latMR,lonMR);
    polesR(i,:)=[lonA,latA];
end

ageMean=poles(:,1);

% polesInput=[1ageMid 2ageSpan 3lonP 4latP 5A95 6Q 7ResNo]
polesInput=[ageMean AgeAll(:,1)-AgeAll(:,2) polesR poles(:,4) poles(:,6) poles(:,5)];

% polesInputNew=[polesInput nj Kj]=[1ageMid 2ageSpan 3lonP 4latP 5A95 6Q 7ResNo 8nj 9Kj];
% function polesInputNew=wAPWP_4KnPangea(polesInput)
polesInputNew=wAPWP_4KnPangea(polesInput);

% function IndexT=wAPWP_AgeWinAgeErr(ageMid,ageSpan,Tinv,window)
% IndexT=[ageRW,counter,counterSTART,counterEND]
[IndexT,IndexTtmp1,IndexTtmp2,IndexTfin]=wAPWP_AgeWinAgeErr(polesInputNew(:,1),polesInputNew(:,2),Tinv,window);

A95c=10; Wabc=ones(3,1); Wscale=2;
% meanW=[1lonM,2latM,3e95a,4e95b,5omega,6Kx,7Ky,8N];
% polesInputNew=[1ageMid 2ageSpan 3lonP 4latP 5A95 6Q 7ResNo 8nj 9Kj];
polesInputNewW=nan(length(polesInputNew(:,1)),1);
for k=1:length(IndexTfin)
    npts=IndexTfin(k).counter;
    AgeWinAgeErrInd=IndexTfin(k).ind;
    
    %------------------------------------------------------------------
    Wj=[];
    for jj=1:length(AgeWinAgeErrInd)
        j=AgeWinAgeErrInd(jj);
        ageRangeIN=[polesInputNew(j,1)-polesInputNew(j,2)/2 ...
            polesInputNew(j,1)+polesInputNew(j,2)/2];
        windowIN=[IndexT(k,1)-window/2 IndexT(k,1)+window/2];
        % function weight=wAPWP_weight(ageRangeIN,windowIN,A95,A95c,Q,Wabc,normFlag)
        weight=wAPWP_weight(ageRangeIN,windowIN,polesInputNew(j,5),A95c,polesInputNew(j,6),Wabc,'AgeA95Q');
        Wj=[Wj; weight.Wj];
    end
    %------------------------------------------------------------------
    polesInputNewW(AgeWinAgeErrInd)=Wj;
    
    polesInputNewIN=polesInputNew(AgeWinAgeErrInd,:);
    
    % change the scale of weight to change the size of error ellipses
    polesInputNewIN(:,6)=Wj*Wscale;
    
    meanW(k)=wAPWP_meanW(polesInputNewIN);
    meanW(k).ageWin=[IndexT(k,1) IndexT(k,1)-window/2 IndexT(k,1)+window/2];
end

polesInputNew=[polesInputNew polesInputNewW];
RMempty=[];
for k=1:length(IndexTfin)
    AgeWinAgeErrInd=IndexTfin(k).ind;
    meanW(k).sumQ=sum(polesInputNew(AgeWinAgeErrInd,6));
    if isempty(meanW(k).ageAll)==1 RMempty=[RMempty;k]; end
end
meanW(RMempty)=[];

%% plot results
ageFilter=[250 500];
ageSelectU=[ageFilter(1):10:ageFilter(2)]';
cmapstep=10;
agesUstep=floor(ageSelectU/cmapstep)*cmapstep; agesUstepU=unique(agesUstep);
cmap=jet(length(agesUstepU));

% extract the 300-460 Ma segment
apwpSeg=[300 460];
AgeCont=cell2mat({meanW.ageWin}'); IndAgeCont=find(AgeCont(:,1)>=min(apwpSeg) & AgeCont(:,1)<=max(apwpSeg)); IndAge=IndAgeCont;

flatRGB=[0 .45 .74; .85 .33 .1; .93 .69 .13];

figctrl=figure;
% [left bottom width height]
figctrl.Position=[79 300 722*1.5 484*1.5];

lat_lim = [-90 90]; lon_lim = [-180 180];

lonshift=3;
marksize=10;
A95transp=.05; 

ax=subplot(1,2,1); lbwh=get(ax, 'position');
set(ax,'position',[lbwh(1)-.11,lbwh(2)-.1,lbwh(3)+.14,lbwh(4)+.14]);
axesm(projtype,'Origin',viewA,...
    'MeridianLabel','off','ParallelLabel','off',...
    'MapLatLimit',lat_lim,'MapLonLimit',lon_lim);
framem on; gridm on; axis off; tightmap
coast=load('coast');
geoshow(coast.lat,coast.long,'DisplayType','polygon',...
    'edgecolor',.7*ones(1,3),'FaceColor','none'); 
% geoshow(coast.lat,coast.long,'DisplayType','polygon',...
%     'FaceColor',.85*ones(1,3),'edgecolor','none'); 

if meanEflag==1
    for j=1:length(IndAge)
        i=IndAge(j);
        [elat1,elon1]=ellipse1(meanW(i).meanW(2),meanW(i).meanW(1),...
            [meanW(i).meanW(3) axes2ecc(meanW(i).meanW(3),meanW(i).meanW(4))], meanW(i).meanW(5));
        geoshow(elat1,elon1,'DisplayType','polygon',...
            'edgecolor',.7*ones(1,3),'FaceColor',[.85 .1 1]); alpha(A95transp);
    end
end

for j=1:length(IndAge)
    i=IndAge(j);
    % meanW=[lonM,latM,e95a,e95b,omega,Kx,Ky,npts];
    geoshow(meanW(i).meanW(2),meanW(i).meanW(1),'DisplayType','point','Marker','S',...
        'MarkerSize',meanW(i).sumQ/4,'MarkerFaceColor',[.85 .7 1],'MarkerEdgeColor','k'); 
    if mod(meanW(i).ageWin(1)/10,2)==1  % odd ten-year
        textm(meanW(i).meanW(2),meanW(i).meanW(1)+lonshift,num2str(meanW(i).ageWin(1)),'color',.5*[1 0 1]);
    end
    % apwp4kin=[Age,lonP,latP,A95,N,Qsum];
    apwp4kin(j,:)=[meanW(i).ageWin(1), meanW(i).meanW(1:2), sqrt((meanW(i).meanW(3))^2+(meanW(i).meanW(4))^2), ...
        meanW(i).meanW(end), meanW(i).sumQ(1)];
    % rm4table=[Age,lonM,latM,e95a,e95b,omega,Kx,Ky,npts];
    rm4table(j,:)=[meanW(i).ageWin(1), meanW(i).meanW];
end
%--------------------------------------------------------------------------
latGCtmpW=[]; lonGCtmpW=[]; latGCtmpF=[]; lonGCtmpF=[];
for j=1:length(IndAge)-1
    i=IndAge(j);
    [latGCtmp1,lonGCtmp1] = track2(meanW(i+1).meanW(2),meanW(i+1).meanW(1),...
        meanW(i).meanW(2),meanW(i).meanW(1),[],'degrees',100);
    [latGCtmp2,lonGCtmp2] = track2(meanW(i+1).fishM(2),meanW(i+1).fishM(1),...
        meanW(i).fishM(2),meanW(i).fishM(1),[],'degrees',100);
    latGCtmpW=[latGCtmpW; flipud(latGCtmp1)];
    lonGCtmpW=[lonGCtmpW; flipud(lonGCtmp1)];
    latGCtmpF=[latGCtmpF; flipud(latGCtmp2)];
    lonGCtmpF=[lonGCtmpF; flipud(lonGCtmp2)];    
end
resPLOTw=[latGCtmpW lonGCtmpW]; resPLOTf=[latGCtmpF lonGCtmpF];
%--------------------------------------------------------------------------
hMw=geoshow(resPLOTw(:,1),resPLOTw(:,2),'DisplayType','line',...
    'linestyle','-','Color',[.75 .0 .75],'LineWidth',2);

age_intp=[apwpSeg(1):1:apwpSeg(end)]';
ind10Myr=find(mod(age_intp,10)==0);

% do spline
% function RMsmoothed=sphsplW(age,lonP,latP,Q,age_intp,S)
% apwp4kin=[Age,lonP,latP,A95,N,Qsum];
RMsmoothed1Myr=sphsplW(apwp4kin(:,1),apwp4kin(:,2),apwp4kin(:,3),apwp4kin(:,4),...
    age_intp,Smooth_Parameter);
RMsmoothed10Myr=RMsmoothed1Myr(ind10Myr,:); apwp4kin=[apwp4kin RMsmoothed10Myr(:,2:3)];
hMw=geoshow(RMsmoothed1Myr(:,3),RMsmoothed1Myr(:,2),'DisplayType','line',...
    'linestyle','-','Color',[0 .5 0],'LineWidth',2);
geoshow(RMsmoothed10Myr(:,3),RMsmoothed10Myr(:,2),'DisplayType','point','Marker','o',...
    'MarkerSize',10,'MarkerFaceColor',[.76 .87 .78],'MarkerEdgeColor','k');
for i=1:length(RMsmoothed10Myr(:,1))
    if mod(RMsmoothed10Myr(i,1),20)==0  % multiples of 20
        textm(RMsmoothed10Myr(i,3),RMsmoothed10Myr(i,2)+lonshift,...
            num2str(RMsmoothed10Myr(i,1)),'color',.3*[0 1 0]);
    end
end

if isnan(apwpSeg)==0
    title('\fontsize{14} APW path')
end

% projtype='mollweid';
ax=subplot(1,2,2); lbwh=get(ax, 'position');
set(ax,'position',[lbwh(1)-.11,lbwh(2)-.1,lbwh(3)+.14,lbwh(4)+.14]);
axesm(projtype,'Origin',viewA,...
    'MeridianLabel','off','ParallelLabel','off',...
    'MapLatLimit',lat_lim,'MapLonLimit',lon_lim);
framem on; gridm on; axis off; tightmap
coast=load('coast');
% geoshow(coast.lat,coast.long,'DisplayType','polygon',...
%     'edgecolor',.85*ones(1,3),'FaceColor','none'); 
geoshow(coast.lat,coast.long,'DisplayType','polygon',...
    'FaceColor',.85*ones(1,3),'edgecolor','none'); 

% plot A95 of paleopoles
for i=1:length(polesInputNew(:,1))
    [minA,indexA] = min(abs(agesUstepU-floor(polesInputNew(i,1)/cmapstep)*cmapstep));
    [elat1,elon1]=scircle1(polesInputNew(i,4),polesInputNew(i,3),polesInputNew(i,5));
    geoshow(elat1,elon1,'DisplayType','polygon',...
        'edgecolor',cmap(indexA,:),'FaceColor','none'); 
end

% plot paleopoles
for i=1:length(polesInputNew(:,1))
    
    % function weight=wAPWP_weight(ageRangeIN,windowIN,A95,A95c,Q,Wabc,normFlag)
    weight4plot=wAPWP_weight([0 1],[0 1],polesInputNew(i,5),A95c,polesInputNew(i,6),Wabc,'A95Q');
    marksize=weight4plot.Wj*15;
    
    [minA,indexA] = min(abs(agesUstepU-floor(polesInputNew(i,1)/cmapstep)*cmapstep));
    geoshow(polesInputNew(i,4),polesInputNew(i,3),'DisplayType','point','Marker','o',...
        'MarkerSize',marksize,'MarkerFaceColor',cmap(indexA,:),'MarkerEdgeColor','k');
end
for i=1:length(polesInputNew(:,1))
    if AgeResNoflag==0
        textm(polesInputNew(i,4),polesInputNew(i,3)+lonshift,num2str(i),'color','k');
    elseif AgeResNoflag==1
        textm(polesInputNew(i,4),polesInputNew(i,3)+lonshift,num2str(polesInputNew(i,1)),'color','k');
    elseif AgeResNoflag==2
        textm(polesInputNew(i,4),polesInputNew(i,3)+lonshift,num2str(polesInputNew(i,7)),'color','k');
    elseif AgeResNoflag==3
    end
end

if isempty(agesUstepU)==0
    caxis([min(agesUstepU) max(agesUstepU)]);
end
colormap(cmap); colorbar;

if isnan(apwpSeg)==0
    title('\fontsize{14} Paleopoles')
end

% output=[1Age 2lonRM 3latRM 4e95a 5e95b 6omega 7Kx 8Ky 9N 10lonSP 11latSP]
output=[rm4table apwp4kin(:,7:8)]
