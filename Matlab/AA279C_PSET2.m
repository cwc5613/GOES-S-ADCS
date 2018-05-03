clear;
clc;
close all;

m = 83;
A_tot = 4.11;
cm = [0;0.39;-0.01];

Ix = [0,0.95,-0.33];
Iy = [0,0.33,0.95];
Iz = [1,0,0];

Px = 5.61;
Py = 7.03;
Pz = 7.49;

J = [7.49 0 0;
    0 5.76 -0.44;
    0 -0.44 6.88];

J_diag = diag(J);
W = [10;0;0];

%% Euler's Equation

%[t,y] = ode45(@(t,y) euler_eqn_solv(J,mu), tspan, initCond, opts);

dt = 0.5;
for t = 1:dt:100
    k_w(:,1) = dt*getWdot(W,t,J_diag);
    k_w(:,2) = dt*getWdot(W+.5*k_w(:,1),t+.5*dt,J_diag);
    k_w(:,3) = dt*getWdot(W+.5*k_w(:,2),t+.5*dt,J_diag);
    k_w(:,4) = dt*getWdot(W+k_w(:,3),t+dt,J_diag);
    W = W+(1/6)*(k_w*[1;2;2;1]);
end

%% Safe Mode



%% Spacecraft Dynamics


%% Helper Functions

function wd = euler_eqn_solv(J,mu)
%wd(1) = [r(4:6); -mu*r(1:3)/norm(r(1:3))^3];
end

function [Wdot] = getWdot( W, T, I )
Imat = diag(I);
Imat_inv = diag(1./I);
Wdot = Imat_inv*cross(Imat*W,W);
end


