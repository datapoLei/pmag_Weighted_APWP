function [IndexT,tmp1,tmp2,IndexTfin]=wAPWP_AgeWinAgeErr(ageMid,ageSpan,Tinv,window)
% wAPWP_AgeWinAgeErr(ageMid,ageSpan,Tinv,window)
% IndexT=[ageRW,counter,counterSTART,counterEND]
% Index=IndexT(:,2);

% N=length(ageMid);
winHalf=window/2;
ageMax=ageMid+ageSpan/2; ageMaxMax=max(ageMax);
ageMin=ageMid-ageSpan/2; ageMinMin=min(ageMin);

%--------------------------------------------------------------------------
% ensure that the age vector starts with a nearest multiple of 10
ageRW=[floor((ageMinMin-winHalf)/winHalf)*winHalf : Tinv : ceil((ageMaxMax+winHalf)/winHalf)*winHalf]';
if mod(floor((ageMinMin-winHalf)/winHalf)*winHalf,10)~=0
    ageRW=ageRW-mod(floor((ageMinMin-winHalf)/winHalf)*winHalf,10);
end
%--------------------------------------------------------------------------

IndexT=NaN(length(ageRW),4);
IndTi=NaN;

%--------------------------------------------------------------------------
% discretize age windows (ageRWrange) and age spans of paleopoles (ageSpan)
% export intermediate values
for j=1:length(ageRW)
    ageRWrange=[ageRW(j)-window/2 : .1 : ageRW(j)+window/2]';
    tmp1(j).ageRWrange=ageRWrange; 
    tmp1(j).ageRWrangeBound=[ageRW(j)-window/2 ageRW(j)+window/2];
end
for i=1:length(ageMid)
    ageRange=[ageMid(i)-ageSpan(i)/2 : .1 : ageMid(i)+ageSpan(i)/2]';
    tmp2(i).ageRange=ageRange;
    tmp2(i).ageRangeBound=[ageMid(i)-ageSpan(i)/2 ageMid(i)+ageSpan(i)/2];
end
%--------------------------------------------------------------------------

for j=1:length(ageRW)
    counter=0;
    tmp3tmp=[];
    for i=1:length(ageMid)      
        if isempty(intersect(tmp1(j).ageRWrange, tmp2(i).ageRange))==0 & ...
                length(intersect(tmp1(j).ageRWrange, tmp2(i).ageRange))>1

            counter=counter+1; IndTi=i;
            tmp3tmp=[tmp3tmp; i];
        end 
    end
    
    Index(j)=counter;
    % IndexT=[ageRW,counter,counterSTART,counterEND]
    IndexT(j,:)=[ageRW(j),counter,IndTi-counter+1,IndTi];
    
    %--------------------------------------------------------------------------
    % for debug: IndexTfin - the final results of age grouping
    IndexTfin(j).ageRW=[ageRW(j) ageRW(j)-window/2 ageRW(j)+window/2]; 
    IndexTfin(j).counter=counter;
    IndexTfin(j).ind=tmp3tmp;
    %--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% for debug: IndexTfin - the final results of age grouping 2019.09.02
negInd=find(ageRW<0);
IndexT(negInd,:)=[]; tmp1(negInd)=[]; tmp2(negInd)=[]; 
IndexTfin(negInd)=[]; 
%--------------------------------------------------------------------------




