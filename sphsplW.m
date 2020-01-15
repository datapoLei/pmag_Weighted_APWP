% do spline
function RMsmoothed=sphsplW(age,lonP,latP,Q,age_intp,S)


% % wrap the data to sphere
% lonP=unwrap(lonP./(180/pi),pi).*(180/pi);
% wrapTo180

% Interpolate data with weight using cubic smoothing spline
lonP_intp=csaps(age,lonP,1/S,age_intp,Q);
latP_intp=csaps(age,latP,1/S,age_intp,Q);

% wrap the data to sphere
% [lonPW_intp,latPW_intp]=toDegrees('radians',lonPW_intp,latPW_intp);
lonPwrap_intp=(lonP_intp);
latPwrap_intp=latP_intp;

% % % % [x,y,z] = sph2cart(lonP/180*pi,latP/180*pi,1);
% % % % 
% % % % % Interpolate data with weight using cubic smoothing spline
% % % % x_intp=csaps(age,x,1/S,age_intp,Q);
% % % % y_intp=csaps(age,y,1/S,age_intp,Q);
% % % % z_intp=csaps(age,z,1/S,age_intp,Q);
% % % % 
% % % % [lonPwrap_intp,latPwrap_intp] = cart2sph(x_intp,y_intp,z_intp);
% % % % 
% % % % RMsmoothed=[age_intp lonPwrap_intp*180*pi latPwrap_intp*180*pi];

RMsmoothed=[age_intp lonPwrap_intp latPwrap_intp];


