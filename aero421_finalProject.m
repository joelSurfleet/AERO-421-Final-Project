% AERO 421 - Project Deliverable 1 
% Group 25

%% Initial Values/Givens
% Mass Properties
clc;clear;

% mass of components
mass_solar = 20; % mass of solar panel, in kg
mass_sensor = 100; % mass of sensor, in kg
mass_bus = 500; % mass of bus, in kg
mass_tot = 640; % satellite stowed mass, in kg

% geometry of components
bus_side = 2; % side length, in meters
solar_height = 2; % side length in x-axis, in meters
solar_base = 3; %  side length in y-axis, in meters
solar_thick = 0.05; % side length in z-axis, in meters
sensor_base = 0.25; % side length in x-axis, in meters
sensor_thick = 0.25; % side length in y-axis, in meters
sensor_height = 1; % side length in z-axis, in meters

% CoM, center of mass;
% assume density of stowed spacecraft is uniform in the 2m cube
CoM_stowed = [0;0;0]; % in m

% deployed CoM, in m
x_deploy = 0;
y_deploy = 0;
z_deploy = mass_sensor*(sensor_height/2+bus_side/2)/mass_tot;
CoM_deploy = [x_deploy;y_deploy;z_deploy];

% J, moment of inertia matrix
% assume the body frame is the principle axis frame 
% stowed moment of inertia, in kg*m^2
Jx_stow = 1/12*mass_tot*(2*bus_side^2); 
Jy_stow = Jx_stow;
Jz_stow = Jx_stow;
J_stow = diag([Jx_stow, Jy_stow, Jz_stow]);

% deployed moment of inertia for sensor, in kg*m^2
I_sensor = diag([1/12*mass_sensor*(sensor_base^2+sensor_height^2); ... Moment of Inertia of Sensor about x axis, m^4
           1/12*mass_sensor*(sensor_base^2+sensor_height^2); ...
           1/12*mass_sensor*(sensor_base^2+sensor_thick^2)]);

rx_sensor = vect2cross([0; 0; bus_side/2 + sensor_height/2 - z_deploy]);

J_sensor = I_sensor - mass_sensor * rx_sensor * rx_sensor;

% deployed moment of inertia for each solar panel, in kg*m^2
I_solar = diag([1/12*mass_solar*(solar_base^2+solar_thick^2);
                1/12*mass_solar*(solar_height^2+solar_thick^2);
                1/12*mass_solar*(solar_height^2+solar_base^2)]);

rx_solar1 = vect2cross([0;bus_side/2+solar_base/2;-z_deploy]);
rx_solar2 = vect2cross([0;-bus_side/2-solar_base/2;-z_deploy]);

J_solar1 = I_solar - mass_solar * rx_solar1 * rx_solar1;
J_solar2 = I_solar - mass_solar * rx_solar2 * rx_solar2;

% deployed moment of inertia for bus, in kg*m^2
I_bus = diag([1/12*mass_bus*(2*bus_side^2); ...
              1/12*mass_bus*(2*bus_side^2); ...
              1/12*mass_bus*(2*bus_side^2)]);

rx_bus = vect2cross([0;0;-z_deploy]);

J_bus = I_bus - mass_bus * rx_bus * rx_bus;

%moment of inertia for deployed configuration, in kg*m^2
J_deploy = J_sensor + J_solar1 + J_solar2 + J_bus

% Orbit Data
mu = 368600; % Km
R_earth = 6378; % Km

COE_0 = [53335.2,0,0,98.43,0,0];

[R_0,V_0] = COE2RV(COE_0,mu);

% Initial Attitude
% Quaternion Relating F_b to F_LVLH
e = [0;0;0];
n = 1;

% Initial Angular Velocity
w_b_ECI_stow   = [-0.05;   0.03;   0.2];
w_b_ECI_deploy = [0.001; -0.001; 0.002];

% Currently, we are going Torque - Free:
T = 0;

% Use Simulink (like ode45) to Create the Plots for the Deliverable

tspan = linspace(0,100,100000);

out = sim("aero421_finalProject_4_10_sim.slx")

%% Part 2
% The spacecraft initial attitude is such that it is aligned with F_LVLH
% Compute initial quaternion and EULER angles relating Fb and Feci

Z_LVLH = -R_0 / norm(R_0);
Y_LVLH = -cross(R_0,V_0) / norm(cross(R_0,V_0));
X_LVLH =  cross(Y_LVLH,Z_LVLH);

C_b_ECI = [X_LVLH;Y_LVLH;Z_LVLH]

n_ECI = (trace(C_b_ECI) + 1) ^ (1/2) / 2
e_ECI = [(C_b_ECI(2,3)-C_b_ECI(3,2)) / (4*n); ...
         (C_b_ECI(3,1)-C_b_ECI(1,3)) / (4*n); ...
         (C_b_ECI(1,2)-C_b_ECI(2,1)) / (4*n)]

phi   = atan(C_b_ECI(2,3) / C_b_ECI(3,2))
theta =-asin(C_b_ECI(1,3))
psi   = atan(C_b_ECI(1,2) / C_b_ECI(1,1))

n_ECI_dot = -1/2 * e_ECI' * w_b_ECI_deploy
e_ECI_dot =  1/2 * (n_ECI * eye(3) + vect2cross(e_ECI)) * w_b_ECI_deploy

eulerRates_ECI = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                  0, cos(phi)           ,-sin(phi);           ...
                  0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
eulerRates_ECI = eulerRates_ECI * w_b_ECI_deploy

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
