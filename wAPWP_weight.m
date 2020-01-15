
function weight=wAPWP_weight(ageRangeIN,windowIN,A95,A95c,Q,Wabc,Wflag)

if nargin<7
    Wflag='AgeA95Q';
end

% Wabc: Weighting coefficients

[W_dt,overlap]=wAPWP_ageOverlap(ageRangeIN,windowIN);

if A95>A95c
    W_A95=A95c/A95;
else
    W_A95=1;
end

W_Q=Q/7;

if strcmpi(Wflag,'AgeA95Q')==1 
    Wabc=Wabc;
elseif strcmpi(Wflag,'Age')==1
    Wabc(2)=0; Wabc(3)=0;
elseif strcmpi(Wflag,'A95')==1
    Wabc(1)=0; Wabc(3)=0;
elseif strcmpi(Wflag,'Q')==1
    Wabc(1)=0; Wabc(2)=0;
elseif strcmpi(Wflag,'AgeA95')==1 
    Wabc(3)=0;
elseif strcmpi(Wflag,'AgeQ')==1 
    Wabc(2)=0;
elseif strcmpi(Wflag,'A95Q')==1 
    Wabc(1)=0;
end

Wj=(Wabc(1)*W_dt+Wabc(2)*W_A95+Wabc(3)*W_Q)/sum(Wabc);

weight.Wj=Wj; 
weight.W_dt=[W_dt,overlap]; weight.W_A95=W_A95; weight.W_Q=W_Q; 
weight.Wabc=Wabc;


function [W_dt,overlap]=wAPWP_ageOverlap(ageRangeIN,windowIN)

% ageRangeIN=[4.2 7.4];
% windowIN=[3.1 12.7];

% ageRangeIN=[4.2 7.4];
% windowIN=[6.1 12.7];
% 
% ageRangeIN=[4.2 5.4];
% windowIN=[6.1 12.7];

ageRange=round([min(ageRangeIN):.1:max(ageRangeIN)],1);
window=round([min(windowIN):.1:max(windowIN)],1);

ages4calc=(intersect(ageRange,window))';
ages4calc=sortrows(ages4calc,'ascend');

if isempty(ages4calc)==0
    overlap=max(ages4calc)-min(ages4calc);
else
    overlap=0;
end

W_dt=overlap/(max(ageRange)-min(ageRange));



