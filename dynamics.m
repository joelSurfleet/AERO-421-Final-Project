function domega = dynamics(T, J, omega)
%#codegen

% Compute the derivative of the angular velocity $\omega$ here

% Using the Euler Equations and T = 0

wXdot = (-(J(3,3)-J(2,2))*omega(2)*omega(3))/J(1,1) + T;
wYdot = (-(J(1,1)-J(3,3))*omega(1)*omega(3))/J(2,2) + T;
wZdot = (-(J(2,2)-J(1,1))*omega(1)*omega(2))/J(3,3) + T;
domega = [wXdot; wYdot; wZdot];
