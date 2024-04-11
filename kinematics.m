function [dq,dE] = kinematics(omega, q, E)
%#codegen

% unpack quaternion

% compute the derivative of the quaternion components here

% compute the derivative of the Euler Angle components here


% The spacecraft initial attitude is such that it is aligned with F_LVLH
% Compute initial quaternion and EULER angles relating Fb and Feci

n = q(1)
e = q(2:4)

dn =-1/2 * e' * omega
de = 1/2 * (n * eye(3) + vect2cross(e)) * omega

dq = [dn; de];

phi   = E(1);
theta = E(2);

dE = [1, sin(phi)*tan(theta), cos(phi)*tan(theta); ...
     0, cos(phi)           ,-sin(phi);            ...
     0, sin(phi)*sec(theta), cos(phi)*sec(theta)] * omega;

end
