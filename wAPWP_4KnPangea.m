
% polesInput=[1ageMid 2ageSpan 3lonP 4latP 5A95 6Q 7ResNo]
% polesInputNew=[polesInput nj Kj]=[1ageMid 2ageSpan 3lonP 4latP 5A95 6Q 7ResNo 8nj 9Kj];

function polesInputNew=wAPWP_4KnPangea(polesInput)

% Assign the same Ns and Ks to all paleopoles
for i=1:length(polesInput(:,1))
    nj(i,1)=10;
    Kj(i,1)=19.6;
end

polesInputNew=[polesInput nj Kj];




