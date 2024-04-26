function [b,m] = magField(rVec)

%% Magnetic Field Model

aMag = 6371.2e3; % m
g0 = -1450.9e-9; % T
h1 = 4652.5e-9; % T
g1 = -29404.8e-9; % T
r = norm(R_0);

m = aMag^3*[g0;h1;g1]; % Dipole Vector For the Magnetic Field m^3*T

% I think that Davey did the COEs (for rECEF)
b = (3*(mehielSat.m'*rVec)*rVec-rVec*mehielSat.m)/r^5;