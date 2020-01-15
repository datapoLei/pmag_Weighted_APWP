
function Ij=wAPWP_pole2I(phi, lambda, K, N)

Rj=wAPWP_pole2R(phi, lambda);

Aj=K*N/(1+K);
Bj=K*N/(1+K);
Cj=2*N/(1+K);

Dj=[Aj 0 0; 0 Bj 0; 0 0 Cj];

Ij=transpose(Rj)*Dj*Rj;


function Rj=wAPWP_pole2R(phi, lambda)

% R11=-sind(lambda)*cosd(phi);
% R12=-sind(lambda)*sind(phi);
% R13=cosd(lambda);

R11=sind(lambda)*cosd(phi);
R12=sind(lambda)*sind(phi);
R13=-cosd(lambda);

R21=-sind(phi);
R22=cosd(phi);
R23=0;

R31=cosd(lambda)*cosd(phi);
R32=cosd(lambda)*sind(phi);
R33=sind(lambda);

Rj=[R11 R12 R13; R21 R22 R23; R31 R32 R33];

