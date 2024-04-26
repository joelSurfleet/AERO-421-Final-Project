% AERO 421 - Project Deliverable 1 
% Group 25

%% Initial Values/Givens
% Mass Properties
clc;clear;

% mass of components
mSen = 100; % mass of sensor, in kg
mSol = 20;  % mass of solar panel, in kg
mBus = 500; % mass of bus, in kg
mTot = 640; % satellite stowed mass, in kg

% dimensions of components
dimSen = [0.25,0.25,1];
dimSol = [2,3,0.05];
dimBus = [2,2,2];

% locations of centroids relative to center of bus
rSen  = [0,0,1.5];
rSol1 = [0,-2.5,0];
rSol2 = [0,2.5,0];
rBus  = [0,0,0];

% CoM, center of mass;
% assume density of stowed spacecraft is uniform in the 2m cube
detumbSat.CoM = [0,0,0]; % in m

% deployed CoM, in m
mehielSat.CoM = [0,0,mSen*(rSen(3))/mTot];

rSen  = rSen  - mehielSat.CoM;
rSol1 = rSol1 - mehielSat.CoM;
rSol2 = rSol2 - mehielSat.CoM;
rBus  = rBus  - mehielSat.CoM;

% Body Geometry:

% mehielSat.names = names of surfaces
% mehielSat.n = Normal Vectors
% mehielSat.c = Centers
% mehielSat.A = Areas

% Done in order of center z coord from -z to +z 
% if 2 surfaces have the same z coord, then it goes from -y to +y
% then -x to +x if both z and y are the same

z = [0,0,-1];
Z = [0,0,1];
y = [0,-1,0];
Y = [0,1,0];
x = [-1,0,0];
X = [1,0,0];

% Nomenclature: Object, Direction
mehielSat.names = ["Bus -z";"-y Panel -z";"+y Panel -z";"-y Panel -y"; ... 
    "-y Panel -x";"-y Panel +x";"Bus -y";"Bus -x";"Bus +x";"Bus +y"; ...
    "+y Panel -x";"+y Panel +x";"+y Panel +y";"-y Panel +z"; ...
    "+y Panel +z";"Bus +z";"Sensor -y";"Sensor -x";"Sensor +x"; ...
    "Sensor -y";"Sensor +z"];

mehielSat.n = [z;z;z;y;x;X;y;x;X;Y;x;X;Y;Z;Z;Z;y;x;X;y;Z];

mehielSat.A = [4;6;6;0.1;0.15;0.15;4;4;4;4;0.15;0.15;0.1;6;6;4;0.25;0.25; ...
    0.25;0.25;0.0625];

mehielSat.C = [[0,0,-1];[0,-2.5,-0.025];[0,2.5,-0.025];[0,-4,0];[-1,-2.5,0]; ...
    [1,-2.5,0];[0,1,0];[-1,0,0];[1,0,0];[0,1,0];[-1,2.5,0];[1,2.5,0]; ...
    [0,4,0];[0,-2.5,0.025];[0,2.5,0.025];[0,0,1];[0,-0.125,1.5]; ...
    [-0.125,0,1.5];[0.125,0,1.5];[0,0.125,1.5];[0,0,2]];

mehielSat.C = mehielSat.C - mehielSat.CoM;

indices = [1,7,8,9,10,16];

detumbSat.n = mehielSat.n(indices,:);
detumbSat.A = mehielSat.A(indices,:);
detumbSat.C = mehielSat.C(indices,:);

% J, moment of inertia matrix
% assume the body frame is the principle axis frame 
% stowed moment of inertia, in kg*m^2
detumbSat.J = J(mTot,dimBus,detumbSat.CoM);

%moment of inertia for deployed configuration, in kg*m^2
mehielSat.J = J(mSen,dimSen,rSen) + J(mSol,dimSol,rSol1) ...
    + J(mSol,dimSol,rSol2) + J(mBus,dimBus,rBus);

% Orbit Data
mu = 398600; % Km

COE_0 = [53335.2,0,0,98.43,0,0];

[R_0,V_0] = COE2RV(COE_0,mu);

Period = 2*pi*sqrt(norm(R_0)^3/mu);

% Initial Attitude
% Quaternion Relating F_b to F_LVLH
% e0 = [0;0;0];
% n0 = 1;

% Initial Angular Velocity
detumbSat.w0 = [-0.05;   0.03;   0.2];
mehielSat.w0 = [0.001; -0.001; 0.002];

%% Part 2
% The spacecraft initial attitude is such that it is aligned with F_LVLH
% Compute initial quaternion and EULER angles relating Fb and Feci

T = 0;

Z_LVLH = -R_0 / norm(R_0);
Y_LVLH = -cross(R_0,V_0) / norm(cross(R_0,V_0));
X_LVLH =  cross(Y_LVLH,Z_LVLH);

C = [X_LVLH;Y_LVLH;Z_LVLH];

n0 = (trace(C) + 1) ^ (1/2) / 2;
e0 = [(C(2,3)-C(3,2)) / (4*n0); ...
      (C(3,1)-C(1,3)) / (4*n0); ...
      (C(1,2)-C(2,1)) / (4*n0)];

q0 = [n0;e0];

phi0   = atan2(C(2,3), C(3,3));
theta0 = -asin(C(1,3));
psi0   = atan2(C(1,2), C(1,1));

E0 = [phi0;theta0;psi0];

%% Magnetic Field Model

aMag = 6371.2e3; % m
g0 = -1450.9e-9; % T
h1 = 4652.5e-9; % T
g1 = -29404.8e-9; % T
rMag = norm(R_0);

mehielSat.m = aMag^3*[g0;h1;g1]; % Dipole Vector For the Magnetic Field m^3*T

% I think that Davey did the COEs (for rECEF)
mehielSat.b = (3*(mehielSat.m'*rMag)*rMag-rMag*mehielSat.m)/rMag^5;

out = sim("aero421_finalProjectSim.slx");

%% Plot Results

out.E_b_ECI(:,2:4) = out.E_b_ECI(:,2:4) .* (180/pi);

close all;

figure('numbertitle','off','name','final project part 2','windowstate','maximized')

sgtitle("Final Project Part 2")

subplot(3,1,1)
grid on; hold on;
plot(out.w_b_ECI(:,1),out.w_b_ECI(:,2))
plot(out.w_b_ECI(:,1),out.w_b_ECI(:,3))
plot(out.w_b_ECI(:,1),out.w_b_ECI(:,4))

legend("\omega_{x}","\omega_{y}","\omega_{z}")
title("Angular Velocities")
xlabel("time (sec)")
ylabel("angular velocity (rad/s)")

subplot(3,1,2)
grid on; hold on;
plot(out.q_b_ECI(:,1),out.q_b_ECI(:,2))
plot(out.q_b_ECI(:,1),out.q_b_ECI(:,3))
plot(out.q_b_ECI(:,1),out.q_b_ECI(:,4))
plot(out.q_b_ECI(:,1),out.q_b_ECI(:,5))

legend("\eta","\epsilon_{1}","\epsilon_{2}","\epsilon_{3}")
title("Quaternions")
xlabel("time (sec)")
ylabel("Quaternion Parameter")

subplot(3,1,3)
grid on; hold on;
plot(out.E_b_ECI(:,1),out.E_b_ECI(:,2))
plot(out.E_b_ECI(:,1),out.E_b_ECI(:,3))
plot(out.E_b_ECI(:,1),out.E_b_ECI(:,4))

legend("\phi","\theta","\psi")
title("Euler Angles")
xlabel("time (sec)")
ylabel("angle (degrees)")

%% Functions

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

function ax = vect2cross(a)

ax = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];

end

function J = J(m,dim,r)
x = dim(1);
y = dim(2);
z = dim(3);

I = diag([m*(y^2+z^2); ...
          m*(x^2+z^2); ...
          m*(x^2+y^2)]);

rx = vect2cross(r);

J = I/12 - m * rx * rx;
end