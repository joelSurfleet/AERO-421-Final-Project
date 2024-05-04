% AERO 421 - MehielSat (Final Project)
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
dimSen = [0.25;0.25;1];
dimSol = [2;3;0.05];
dimBus = [2;2;2];

% locations of centroids relative to center of bus
rSen  = [0;0;1.5];
rSol1 = [0;-2.5;0];
rSol2 = [0;2.5;0];
rBus  = [0;0;0];

% CoM, center of mass;
% assume density of stowed spacecraft is uniform in the 2m cube
detumbSat.CoM = [0;0;0]; % in m

% deployed CoM, in m
mehielSat.CoM = [0;0;mSen*(rSen(3))/mTot];

% distance between mehielSat CoM and Component center
rSen  = rSen  - mehielSat.CoM;
rSol1 = rSol1 - mehielSat.CoM;
rSol2 = rSol2 - mehielSat.CoM;
rBus  = rBus  - mehielSat.CoM;

% Body Geometry:

% Nomenclature:
% mehielSat.names = names of surfaces
% mehielSat.n = Normal Vectors
% mehielSat.c = Centers
% mehielSat.A = Areas

% Methodology
% Done in order of center z coord from -z to +z 
% if 2 surfaces have the same z coord, then it goes from -y to +y
% then -x to +x if both z and y are the same

% Unit vectors just to make things easy
z = [0;0;-1];
Z = [0;0;1];
y = [0;-1;0];
Y = [0;1;0];
x = [-1;0;0];
X = [1;0;0];

% List of Surface names in order
% Nomenclature: Object, Direction
mehielSat.names = ["Bus -z","-y Panel -z","+y Panel -z","-y Panel -y", ... 
    "-y Panel -x","-y Panel +x","Bus -y","Bus -x","Bus +x","Bus +y", ...
    "+y Panel -x","+y Panel +x","+y Panel +y","-y Panel +z", ...
    "+y Panel +z","Bus +z","Sensor -y","Sensor -x","Sensor +x", ...
    "Sensor -y","Sensor +z"];

% Normal Vectors of mehielsat surfaces
mehielSat.n = [z,z,z,y,x,X,y,x,X,Y,x,X,Y,Z,Z,Z,y,x,X,y,Z];

% Areas of surfaces
mehielSat.A = [4,6,6,0.1,0.15,0.15,4,4,4,4,0.15,0.15,0.1,6,6,4,0.25,0.25, ...
    0.25,0.25,0.0625];

% Centers of surfaces relative to center of bus
mehielSat.C = [[0;0;-1],[0;-2.5;-0.025],[0;2.5;-0.025],[0;-4;0],[-1;-2.5;0], ...
    [1;-2.5;0],[0;1;0],[-1;0;0],[1;0;0],[0;1;0],[-1;2.5;0],[1;2.5;0], ...
    [0;4;0],[0;-2.5;0.025],[0;2.5;0.025],[0;0;1],[0;-0.125;1.5], ...
    [-0.125;0;1.5],[0.125;0;1.5],[0;0.125;1.5],[0;0;2]];

% Adjust surface centers for COM
mehielSat.C = mehielSat.C - mehielSat.CoM;

% Bus inidices
indices = [1,7,8,9,10,16];

detumbSat.n = mehielSat.n(:,indices);
detumbSat.A = mehielSat.A(:,indices);
detumbSat.C = mehielSat.C(:,indices);

% J, moment of inertia matrix
% assume the body frame is the principle axis frame 
% stowed moment of inertia, in kg*m^2
detumbSat.J = J(mTot,dimBus,detumbSat.CoM);

% Moment of Inertia for deployed configuration, in kg*m^2
mehielSat.J = J(mSen,dimSen,rSen) + J(mSol,dimSol,rSol1) ...
    + J(mSol,dimSol,rSol2) + J(mBus,dimBus,rBus);

% Orbital Data

mu = 398600; % gravitational parameter, in km
COE_0 = [53335.2,0,0,98.43,0,0]; % initial orbital elements [h,ecc,RAAN,inc,w,omega]
[R_0,V_0] = COE2RV(COE_0,mu); % initial R and V vectors in LVLH
state_0 = [R_0';V_0']; % inital state 
Period = 2*pi*sqrt(norm(R_0)^3/mu); % period of mehielSat orbit

% Initial Angular Velocity
detumbSat.w0 = [-0.05;   0.03;   0.2];
mehielSat.w0 = [0.001; -0.001; 0.002];

% Magnetic Dipole
mehielSat.m = [0;0;-0.5]; % A*m^2

% Part 2
% The spacecraft initial attitude is such that it is aligned with F_LVLH
% Compute initial quaternion and EULER angles relating Fb and Feci
Z_LVLH = -R_0 / norm(R_0);
Y_LVLH = -cross(R_0,V_0) / norm(cross(R_0,V_0));
X_LVLH =  cross(Y_LVLH,Z_LVLH);

C = [X_LVLH;Y_LVLH;Z_LVLH];

omega_LVLH0 = mehielSat.w0 + 2*pi/Period*Y_LVLH';

% Initial Attitude of MehielSat
% inital Quaternion Relating F_b to F_LVLH
q0_LVLH = [0;0;0;1];

% initial Euler Angles relating F_b to F_LVLH
E0_LVLH = [0;0;0];

n0 = (trace(C) + 1) ^ (1/2) / 2;
e0 = [(C(2,3)-C(3,2)) / (4*n0); ...
      (C(3,1)-C(1,3)) / (4*n0); ...
      (C(1,2)-C(2,1)) / (4*n0)];

q0 = [e0;n0];

phi0   = atan2(C(2,3), C(3,3));
theta0 = -asin(C(1,3));
psi0   = atan2(C(1,2), C(1,1));

E0 = [phi0;theta0;psi0];

s0 = X; % ECI
%% Run Sim
out = sim("Fries.slx");

%% Plot Results

out.E(:,2:4) = out.E(:,2:4) .* (180/pi);
out.E_LVLH(:,2:4) = out.E_LVLH(:,2:4) .* (180/pi);

close all;

figure('numbertitle','off','name','final project part 4','windowstate','maximized')

sgtitle("Spacecraft Attitude over 5 Periods")

subplot(3,1,1)
grid on; hold on;
plot(out.w(:,1),out.w(:,2))
plot(out.w(:,1),out.w(:,3))
plot(out.w(:,1),out.w(:,4))

legend("\omega_{x}","\omega_{y}","\omega_{z}")
title("Angular Velocities")
xlabel("time (sec)")
ylabel("angular velocity (rad/s)")

subplot(3,1,2)
grid on; hold on;
plot(out.q(:,1),out.q(:,2))
plot(out.q(:,1),out.q(:,3))
plot(out.q(:,1),out.q(:,4))
plot(out.q(:,1),out.q(:,5))

legend("\eta","\epsilon_{1}","\epsilon_{2}","\epsilon_{3}")
title("Quaternions")
xlabel("time (sec)")
ylabel("Quaternion Parameter")

subplot(3,1,3)
grid on; hold on;
plot(out.E(:,1),out.E(:,2))
plot(out.E(:,1),out.E(:,3))
plot(out.E(:,1),out.E(:,4))

legend("\phi","\theta","\psi")
title("Euler Angles")
xlabel("time (sec)")
ylabel("angle (degrees)")

figure('numbertitle','off','name','final project part 4','windowstate','maximized')

sgtitle("Disturbance Torques on Spacecraft over 5 Periods")

subplot(4,1,1)
grid on; hold on;
plot(out.tout(:,1),out.Ta(1,:))
plot(out.tout(:,1),out.Ta(2,:))
plot(out.tout(:,1),out.Ta(3,:))

legend("T_{ax}","T_{ay}","T_{az}")
title("Atmospheric Drag Torque")
xlabel("time (sec)")
ylabel("Torque (N*m)")

subplot(4,1,2)
grid on; hold on;
plot(out.tout(:,1),out.Tb(1,:))
plot(out.tout(:,1),out.Tb(2,:))
plot(out.tout(:,1),out.Tb(3,:))

legend("T_{bx}","T_{by}","T_{bz}")
title("Magnetic Torque")
xlabel("time (sec)")
ylabel("Torque (N*m)")

subplot(4,1,3)
grid on; hold on;
plot(out.tout(:,1),out.Ts(1,:))
plot(out.tout(:,1),out.Ts(2,:))
plot(out.tout(:,1),out.Ts(3,:))

legend("T_{sx}","T_{sy}","T_{sz}")
title("SRP Torque")
xlabel("time (sec)")
ylabel("Torque (N*m)")

subplot(4,1,4)
grid on; hold on;
plot(out.tout(:,1),out.Tg(1,:))
plot(out.tout(:,1),out.Tg(2,:))
plot(out.tout(:,1),out.Tg(3,:))

legend("T_{gx}","T_{gy}","T_{gz}")
title("Gravity Gradient Torque")
xlabel("time (sec)")
ylabel("Torque (N*m)")

figure('numbertitle','off','name','final project part 4','windowstate','maximized')

sgtitle("Spacecraft Attitude in LVLH over 5 Periods")

subplot(3,1,1)
grid on; hold on;
plot(out.tout(:,1),out.w_LVLH(:,2))
plot(out.tout(:,1),out.w_LVLH(:,3))
plot(out.tout(:,1),out.w_LVLH(:,4))

legend("\omega_{x}","\omega_{y}","\omega_{z}")
title("Angular Velocities")
xlabel("time (sec)")
ylabel("angular velocity (rad/s)")

subplot(3,1,2)
grid on; hold on;
plot(out.tout(:,1),out.q_LVLH(:,2))
plot(out.tout(:,1),out.q_LVLH(:,3))
plot(out.tout(:,1),out.q_LVLH(:,4))
plot(out.tout(:,1),out.q_LVLH(:,5))

legend("\eta","\epsilon_{1}","\epsilon_{2}","\epsilon_{3}")
title("Quaternions")
xlabel("time (sec)")
ylabel("Quaternion Parameter")

subplot(3,1,3)
grid on; hold on;
plot(out.tout(:,1),out.E_LVLH(:,2))
plot(out.tout(:,1),out.E_LVLH(:,3))
plot(out.tout(:,1),out.E_LVLH(:,4))

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
%   mu: gravitational parameter of object being orbited

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