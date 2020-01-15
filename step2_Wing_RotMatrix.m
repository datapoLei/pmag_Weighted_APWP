
function RotMat=step2_Wing_RotMatrix(lonE,latE,omega)

% convert to cartesian coordinates:
[Ex,Ey,Ez]=sph2cart(lonE*pi/180,latE*pi/180,1);

% set up rotation matrix R
R11=Ex*Ex*(1-cosd(omega))+cosd(omega);
R12=Ex*Ey*(1-cosd(omega))-Ez*sind(omega);
R13=Ex*Ez*(1-cosd(omega))+Ey*sind(omega);
R21=Ey*Ex*(1-cosd(omega))+Ez*sind(omega);
R22=Ey*Ey*(1-cosd(omega))+cosd(omega);
R23=Ey*Ez*(1-cosd(omega))-Ex*sind(omega);
R31=Ez*Ex*(1-cosd(omega))-Ey*sind(omega);
R32=Ez*Ey*(1-cosd(omega))+Ex*sind(omega);
R33=Ez*Ez*(1-cosd(omega))+cosd(omega);

% combine results
RotMat=[R11 R12 R13;R21 R22 R23;R31 R32 R33];