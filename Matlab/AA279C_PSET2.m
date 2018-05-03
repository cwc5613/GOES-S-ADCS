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

w_max = [10;0;0];
h = floor(sum(J*w_max));
[X,Y,Z] = sphere(100);

%% Euler's Equation

dt = 0.1;
timeStop = 100;
tspan = 0:dt:timeStop;

w01 = [10;0.01;0.01];
w02 = [0.01;10;0.01];
w03 = [0.01;0.01;10];

[t,w] = ode45(@(t,w) euler_eqn_solv(t,w,J), tspan, w03);

figure
hold on
 surf(X*h, Y*h, Z*h);
plot3(w(:,1),w(:,2),w(:,3),'LineWidth',2)
hold off
axis equal

%% Safe Mode



%% Spacecraft Dynamics


%% Helper Functions

function wdot = euler_eqn_solv(t,w,J)

wdot(1) = (-(J(3,3)-J(2,2))*w(2)*w(3))/(J(1,1));
wdot(2) = (-(J(1,1)-J(3,3))*w(1)*w(3))/(J(2,2));
wdot(3) = (-(J(2,2)-J(1,1))*w(2)*w(1))/(J(3,3));

wdot = wdot';

end


