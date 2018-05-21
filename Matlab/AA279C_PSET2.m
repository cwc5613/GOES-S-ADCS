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
    0 5.61 0;
    0 0 7.03];

J_diag = diag(J);

w_max = [10;0;0];
h_max = sum(J*w_max);
h = sum(J*w_max);
[X,Y,Z] = sphere(50);

%% Euler's Equation

dt = 0.1;
timeStop = 15;
tspan = 0:dt:timeStop;

count = 1;
c = 6;

for i = 0:c
    for j = 1:c
        
        w1 = [0;i;j];
        w1(1) = sqrt((10*J(1,1))^2-(J(2,2)*w1(2))^2-(J(3,3)*w1(3))^2)/J(1,1);
        w2 = [i;0;j];
        w2(2) = sqrt((10*J(1,1))^2-(J(1,1)*w2(1))^2-(J(3,3)*w2(3))^2)/J(2,2);
        w3 = [0;i;j];
        w3(3) = sqrt((10*J(1,1))^2-(J(2,2)*w3(2))^2-(J(1,1)*w3(1))^2)/J(3,3);
        
        [t1,wt1] = ode45(@(t,w) euler_eqn_solv(t,w,J), tspan, w1);
        [t2,wt2] = ode45(@(t,w) euler_eqn_solv(t,w,J), tspan, w2);
        [t3,wt3] = ode45(@(t,w) euler_eqn_solv(t,w,J), tspan, w3);
        
        ht1(:,:,count) = J*wt1';
        ht2(:,:,count) = J*wt2';
        ht3(:,:,count) = J*wt3';
        
        wt1(:,:,count) = wt1;
        wt2(:,:,count) = wt2;
        wt3(:,:,count) = wt3;
        
        count = count + 1;
    end
end

weq01 = [0;0;0];
weq01(1) = sqrt((10*J(1,1))^2-(J(2,2)*weq01(2))^2-(J(3,3)*weq01(3))^2)/J(1,1);
weq02 = [0;0;0];
weq02(2) = sqrt((10*J(1,1))^2-(J(1,1)*weq02(1))^2-(J(3,3)*weq02(3))^2)/J(2,2);
weq03 = [0;0;0];
weq03(3) = sqrt((10*J(1,1))^2-(J(2,2)*weq03(2))^2-(J(1,1)*weq03(1))^2)/J(3,3);

heq1 = J*weq01;
heq2 = J*weq02;
heq3 = J*weq03;

CO(:,:,1) = 0.85.*ones(size(X,1)); % red
CO(:,:,2) = 0.85.*ones(size(X,1)); % green
CO(:,:,3) = 0.85.*ones(size(X,1)); % blue

figure
hold on
surf(X*0.99*h, Y*0.99*h, Z*0.99*h,CO,'FaceAlpha',0.95);
lighting phong
shading interp
for i = 1:count-1
    plot3(ht1(1,:,i),ht1(2,:,i),ht1(3,:,i),'r','LineWidth',1)
    plot3(ht2(1,:,i),ht2(2,:,i),ht2(3,:,i),'b','LineWidth',1)
    plot3(ht3(1,:,i),ht3(2,:,i),ht3(3,:,i),'k','LineWidth',1)
    plot3(-ht1(1,:,i),-ht1(2,:,i),-ht1(3,:,i),'r','LineWidth',1)
    plot3(-ht2(1,:,i),-ht2(2,:,i),-ht2(3,:,i),'b','LineWidth',1)
    plot3(-ht3(1,:,i),-ht3(2,:,i),-ht3(3,:,i),'k','LineWidth',1)
end
plot3(heq1(1,:),heq1(2,:),heq1(3,:),'k*','LineWidth',2)
plot3(heq2(1,:),heq2(2,:),heq2(3,:),'k*','LineWidth',2)
plot3(heq3(1,:),heq3(2,:),heq3(3,:),'k*','LineWidth',2)
plot3(-heq1(1,:),-heq1(2,:),-heq1(3,:),'k*','LineWidth',2)
plot3(-heq2(1,:),-heq2(2,:),-heq2(3,:),'k*','LineWidth',2)
plot3(-heq3(1,:),-heq3(2,:),-heq3(3,:),'k*','LineWidth',2)
title(['Angular Momentum Sphere (h = ' num2str(norm(h)) ')'])
axis equal
hold off

%% Safe Mode

%Spin about the negative third axis
w3 = weq03;

%Superspin and Dynamic Balance Rotor Momentum
Jeff = 1.2*J(1,1);
rho = (Jeff-J(3,3))*w3;
rho_s = dot(rho,w3)/norm(w3);
w_s = norm(w3);
w_hat = hat(w3);
rho_test = ([w3 -w_hat]*[w3';w_hat])\...
    [w3 -w_hat]*[w_s*rho_s;cross(-w3,J*w3)];

h_rho = norm(J*w3+rho);

%Superspin and Dynamic Balance Rotor Momentum
[t3,wspin3] = ode45(@(t,w) euler_rotor_solv(t,w,J,rho), tspan, w3);
hspin3 = J*wspin3' + rho;

count = 1;

for i = 0:c
    for j = 1:c
        wdist03 = [i;j;0];
        wdist03(3) = sqrt((h_rho)^2-(J(2,2)*wdist03(2))^2-(J(1,1)*wdist03(1))^2)/J(3,3);
        [t3,wspindist3] = ode45(@(t,w) euler_rotor_solv(t,w,J,rho), tspan, wdist03);
        hspindist3(:,:,count) = J*wspindist3';
        wspindist3(:,:,count) = wspindist3;
        count = count + 1;
    end
end

figure
hold on
surf(X*0.99*h_rho, Y*0.99*h_rho, Z*0.99*h_rho,CO,'FaceAlpha',0.95);
lighting phong
shading interp
plot3(hspin3(1,:),hspin3(2,:),hspin3(3,:),'b*','LineWidth',2)
for i = 1:count-1
    plot3(hspindist3(1,:),hspindist3(2,:),hspindist3(3,:),'k','LineWidth',1)
end
title(['Angular Momentum Sphere (h = ' num2str(norm(h_rho)) ')'])
axis equal
hold off

%% Spacecraft Dynamics

%AA279C_PSET2_orbit;

%% Helper Functions

function wdot = euler_eqn_solv(t,w,J)
wdot = J\-(cross(w,(J*w)));
end

function wdot = euler_rotor_solv(t,w,J,rho)
wdot = J\-(cross(w,(J*w+rho)));
end

function v = unhat(M)
v(1) = -M(2,3);
v(2) = M(1,3);
v(3) = -M(1,2);
end

function M = hat(v)
M = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
end

function x_dot = HSTdynamics(mu,J,rho,init_state)

% Define State x = [r; r_dot; om; quat], x_dot = [r_dot; r_ddot; om_dot; quat_dot]
rvec = init_state(1:3);
vvec = init_state(4:6);
om0 = init_state(7:9);
q0 = init_state(10:13);

v = q0(1:3);
s = q0(4);
qhat = [s*eye(3)+hat(v) v;-v' s];
r = norm(rvec);

% Matrix Linear Equation
x_dot = [zeros(3,3) eye(3,3);
    -(mu/(r^3))*eye(3,3) zeros(3,3)]*[rvec; vvec];

x_dot(7:9,1) = -J\cross(om0,J*om0 + rho);
x_dot(10:13,1) = (1/2)*qhat*[om0; 0];
end

function q = Q2quat(Q)
phi = unhat(logm(Q));
theta = norm(phi);
r = (phi/theta)';
q = [r*sin(theta/2); cos(theta/2)];
end