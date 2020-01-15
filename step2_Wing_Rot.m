
function [lonMR,latMR]=step2_Wing_Rot(lonM,latM,lonE,latE,omega)

% check distance between point M and Euler pole
dist=distance(latM,lonM,latE,lonE);
if min(dist)>90
    [latE,lonE]=antipode(latE,lonE);
    omega=-omega;
end

% convert Point M(lonM,latM) to cartesian coordinates:
[Mx,My,Mz]=sph2cart(lonM*pi/180,latM*pi/180,1);

% calculate the rotation matrix from Euler parameters (lonE,latE,omega)
RotMat=step2_Wing_RotMatrix(lonE,latE,omega);

for i=1:length(lonM)
    % rotate the Point M to new position M'(lonMR,latMR)
    MR(:,i)=RotMat*[Mx(i);My(i);Mz(i)];
end

% coordinates converted back
[lonMR,latMR,rMR]=cart2sph(MR(1,:),MR(2,:),MR(3,:));
lonMR=lonMR*180/pi;
latMR=latMR*180/pi;

% if (lonMR<0)
%     lonMR=lonMR+360;
% end




