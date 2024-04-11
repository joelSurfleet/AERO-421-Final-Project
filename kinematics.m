function [dq,dE] = kinematics(omega, q, E)
%#codegen

% unpack quaternion

% compute the derivative of the quaternion components here

% compute the derivative of the Euler Angle components here


% The spacecraft initial attitude is such that it is aligned with F_LVLH
% Compute initial quaternion and EULER angles relating Fb and Feci

% Orbit Data
mu = 368600; % Km
R_earth = 6378; % Km

COE_0 = [53335.2,0,0,98.43,0,0];

[R_0,V_0] = COE2RV(COE_0,mu);

Z_LVLH = - R_0 / norm(R_0);
Y_LVLH = - cross(R_0,V_0) / norm(cross(R_0,V_0));
X_LVLH = cross(Y_LVLH,Z_LVLH);

C_b_ECI = [X_LVLH;Y_LVLH;Z_LVLH]

n_ECI = (trace(C_b_ECI) + 1) ^ (1/2) / 2
e_ECI = [(C_b_ECI(2,3)-C_b_ECI(3,2)) / (4*n); ...
         (C_b_ECI(3,1)-C_b_ECI(1,3)) / (4*n); ...
         (C_b_ECI(1,2)-C_b_ECI(2,1)) / (4*n)]

phi   = atan(C_b_ECI(2,3) / C_b_ECI(3,2))
theta =-asin(C_b_ECI(1,3))
psi   = atan(C_b_ECI(1,2) / C_b_ECI(1,1))

n_ECI_dot = -1/2 * e_ECI' * w_b_ECI_deploy
e_ECI_dot = 1/2 * (n_ECI * eye(3) + vect2cross(e_ECI)) * w_b_ECI_deploy

dq = [n_ECI_dot; e_ECI_dot];

eulerRates_ECI = [1, sin(phi)*tan(theta), cos(phi)*tan(theta); ...
                  0, cos(phi)           ,-sin(phi);            ...
                  0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
dE = eulerRates_ECI * w_b_ECI_deploy;

end

function [R,V] = COE2RV(coe,mu)

% Converts the 6 classical orbital elements into a position and velocity
% vector
%
% input:
%   coe: an array with the classical orbital elemnts:
%     (1) h     = Angular Velocity
%     (2) ecc   = Eccentricity
%     (3) RAAN  = Right Ascension of ascending node
%     (4) inc   = Inclination in degrees
%     (5) w     = Argument of Perigee in degrees
%     (6) Theta  = True Anomaly in degrees
%   mu: gravitational parameter of object being orbitted

h     = coe(1);
ecc   = coe(2);
RAAN  = coe(3);
inc   = coe(4);
w     = coe(5);
theta = coe(6);

% Find r and v in the perifocal frame
R_pf = (h^2/mu) * (1/(1+ecc*cosd(theta))) * (cosd(theta)*[1;0;0]+sind(theta)*[0;1;0]);
V_pf = (mu/h) * (-sind(theta)*[1;0;0] + (ecc+cosd(theta))*[0;1;0]);

% Perfom euler sequence to get the rotation matrix from perifocal to ECI
C = (basicrot(3,w)*basicrot(1,inc)*basicrot(3,RAAN))';

R = C*R_pf;
R = R';

V = C*V_pf;
V = V';

end

function [C] = basicrot(a,theta)

if a == 1
    C = [1 0 0 ; 0 cosd(theta) sind(theta); 0 -sind(theta) cosd(theta)];
elseif a == 2
    C = [cosd(theta) 0 -sind(theta); 0 1 0; sind(theta) 0 cosd(theta)];
elseif a == 3
    C = [cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1];
end
end