function T = T(sat,m,b,V,rb,J)

% Atmospherice Drag Torgue
Ta = 0;

v = norm(V);

vhat = V./v;

D = 0; % ???????

for i = 1:16
    if dot(sat.n(i,:),vhat) < 0
        continue
    end

    Cp = sat.c(i,:);
    Fa = 1/2 * D * v^2 * 2.5 * sat.A(i,:) * vhat;
    Ta = Ta + cross(Cp,Fa);
end

% Magnetic Field Torque
Tb = cross(m,b);

% Solar Pressure Torque

p = 4.5e-6;

for i = 1:16
    if dot(sat.n(i,:),s) < 0
        continue
    end

    Cp = sat.C(i,:);
    Fa = -p * s * sat.A(i,:);
    Ts = Ts + cross(Cp,Fa);
end

% Gravity Gradient Torque
rx = vect2cross(rb);
Tg = ((3 * mu) / norm(rb)^5) * rx * J * rb;

% Total Torque
T = Ta + Tb + Ts + Tg;

T = 0;

function ax = vect2cross(a)

ax = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];

end

end
