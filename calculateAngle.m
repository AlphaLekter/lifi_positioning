
function angle_deg = calculateAngle(u, v, alpha, beta)
% PARAMETERS
%alpha = roll angle
%beta = yaw angle

d = calculateDistance(u,v);

x = (((u(1)-v(1))/(d))*sin(alpha)*cos(beta)) + (((u(2)-v(2))/(d))*sin(alpha)*sin(beta)) + (((u(3)-v(3))/(d))*cos(beta));

angle = acos(x);

angle_deg = rad2deg(angle);
end