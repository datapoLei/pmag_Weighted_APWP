function [lonE,latE,omega]=step2_Wing_RotM2Euler(RotMat)

% Euler pole longitude
lonE=atan2d(RotMat(1,3)-RotMat(3,1),RotMat(3,2)-RotMat(2,3));
% if (lonE<0)
%     lonE=lonE+360;
% end

% Euler pole latitude
latE=asind((RotMat(2,1)-RotMat(1,2))/...
    (sqrt((RotMat(3,2)-RotMat(2,3))^2+(RotMat(1,3)-RotMat(3,1))^2+...
    (RotMat(2,1)-RotMat(1,2))^2)));

% Finite rotation angle
omega=atan2d(sqrt((RotMat(3,2)-RotMat(2,3))^2+(RotMat(1,3)-RotMat(3,1))^2+...
    (RotMat(2,1)-RotMat(1,2))^2),RotMat(1,1)+RotMat(2,2)+RotMat(3,3)-1);
