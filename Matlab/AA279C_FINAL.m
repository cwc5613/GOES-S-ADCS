%% AA279C Cumulative Final Report
% Author: Christopher Covert
% Course: AA279C

close all,clear all,clc

%% Initialize Various Constant Parameters
% Orbital Elements and Parameters
[a, e, inc, RAAN, w, M0, E, anom, revs_per_day] =...
    TLE_Reader('SKYSAT1_TLE.txt');
mu = 3.986e5; % Earth Standard Gravitational Parameter (km^3/s^2)
T = 2*pi*sqrt((a^3)/mu); % Orbital Period (s)
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
save('Parameters\orbital','a','e','inc','RAAN','w','anom','mu','T','J2','Re')

% Mass Properties, Principle Axes of Inertia
J11 = 7.49;
J22 = 5.61;
J33 = 7.03;
J = diag([J11 J22 J33]); % [x-axis z-axis y-axis]
save('Parameters\mass','J11','J22','J33','J')

% Sensor Specifications
% Star Tracker Specifications
error_ST = 300; % Arcseconds
error_ST = error_ST*(1/3600)*(pi/180); % Radians
V_st = diag((error_ST^2)*ones(1,3));
% Gyro Information Specifications
gyro_rnd = 0.1*(pi/180); % Rate Noise Density, converted from [deg/s/sqrt(Hz)]
gyro_arw = 0.1*(pi/180); % Angle Random Walk, converted from [deg/sqrt(hr)]
W_rnd = (gyro_rnd^2)*eye(3); % Gyro Covariance
W_arw = (gyro_arw^2)*eye(3); % Bias Covariance
W_gyro = blkdiag(W_rnd, W_arw);
% Magnetometer Specifications
error_mag = 4*(pi/180); % Radians
V_mag = (error_mag^2)*eye(3);
save('Parameters\sensor','error_ST','V_st','gyro_rnd','gyro_arw','W_gyro','error_mag','V_mag')

% Compute Environmental Disturbance Parameters
n1 = [1 0 0]'; n2 = [-1 0 0]'; n3 = [0 1 0]'; n4 = [0 -1 0]'; n5 = [0 0 1]'; n6 = [0 0 -1]';
n = [n1 n2 n3 n4 n5 n6]; % Surface Normal Vectors
S = [0.58,0.58,0.96,0.72,0.57,0.58]; %m^2
c = [-0.25,0.07,-0.03;...
    -0.25,0.07,-0.03;...
    0,0.41,-0.36;...
    0,0.01,-0.29;...
    0,0.07,0.22;...
    0,0.08,-0.28]'; %m n = [n1 n2 n3 n4 n5 n6]; % Surface Normal Vectors
params = [S;n;c];

% Reaction Wheel Actuator Jacobian
Bw = [-1/sqrt(3),  1/sqrt(3),  1/sqrt(3), -1/sqrt(3);
       1/sqrt(3),  1/sqrt(3),  1/sqrt(3),  1/sqrt(3);
      -1/sqrt(3), -1/sqrt(3),  1/sqrt(3),  1/sqrt(3)];
save('Parameters\geometry','S','n','c','params','Bw')

%% Superspin and Dynamic Balance (All Calcs Done in Principal Frame)
clc,clearvars,close all
load('Parameters\orbital')
load('Parameters\mass')

w_max = [10;0;0];
h_max = sum(J*w_max);
h = sum(J*w_max);
[X,Y,Z] = sphere(50);

% Euler's Equation

dt = 0.1;
timeStop = 15;
tspan = 0:dt:timeStop;

count = 1;
c = 4;

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

% Safe Mode

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

%% Orbital Dynamics Simulation
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')

tol = 1e-6;

[r_eci, v_eci] = OE2ECI(a, e, deg2rad(inc), deg2rad(RAAN), deg2rad(w), anom, mu);

initCond = [r_eci; v_eci];

orbitCount = 500;
stopTime = orbitCount*T;
stepSize = 100;
tspan = [0:stepSize:stopTime];
opts  = odeset('reltol', tol, 'abstol', tol);

[t,y] = ode45(@(t,y) orb_prop(y,mu), tspan, initCond, opts);

figure
hold on
plot3(y(:,1),y(:,2),y(:,3),'c','LineWidth',1)
earthPlot(1)
axis equal
hold off
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

[t2,y2] = ode45(@(t,y) J2_prop(y,mu,J2,Re), tspan, initCond, opts);

figure
hold on
plot3(y2(:,1),y2(:,2),y2(:,3),'c','LineWidth',1)
earthPlot(1)
axis equal
hold off
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

%% Attitude Dynamics Simulation
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')

% From Superspin
rho = [0;0;20.8612];

% Initial Attitude Conditions
q0 = [0;0;0;1]; % Identity Quaternion
om0 = [0;1;10.7]; % rad/s
x0 = [om0;q0];

% ODE45 Solver for Attitude Dynamics
dt = 0.05;
tspan = 0:dt:15;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)SkySat_Dynamics_Sim(t,x,mu,J,rho,0), tspan, x0, opts);
om_true = x(:,1:3); q_true = x(:,4:7);

% Plot Quaternion
figure
subplot(4,1,1),plot(t,q_true(:,1)),ylabel('q_1'),title('SkySat1 Quaternion Dynamics')
subplot(4,1,2),plot(t,q_true(:,2)),ylabel('q_2')
subplot(4,1,3),plot(t,q_true(:,3)),ylabel('q_3')
subplot(4,1,4),plot(t,q_true(:,4)),ylabel('q_4'),xlabel('Time [sec]')

% Rerun simulation for longer duration for pointing error
dt = 0.05;
tspan = 0:dt:60;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)SkySat_Dynamics_Sim(t,x,mu,J,rho,0), tspan, x0, opts);
om_true = x(:,1:3); q_true = x(:,4:7);

for i = 1:size(q_true,1)
    err(:,i) = abs(quat2phi(q_true(i,:)))-[0;0;1];
end

% Plot pointing error
figure
plot(t,vecnorm(err));...
    ylabel('Pointing Error (deg)'),xlabel('Time [sec]'),...
    title('SkySat1 Pointing Error'),...
    ylim([0,10]);

%% Static Attitude Estimation (Using Magnetometer + GPS)
load('Parameters\sensor')

% Noisy Magnetometer Body-Frame Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = generateNoise(x,V_mag,1,rN);

% Perform Static Estimation at Each Time Step
for ii = 1:length(x)
    rB = [rB_noisy(1:3,ii) rB_noisy(4:6,ii)];   
    [q_triad(:,ii),~] = triad(rN,rB);
    [q_dav(:,ii),~] = qdaven(rN,rB,error_mag);
    [q_svd(:,ii),~] = qsvd(rN,rB,error_mag);
    q_triad(:,ii) = sign(x(ii,4:7)').*abs(q_triad(:,ii));
    error_triad(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_triad(:,ii)));
    q_dav(:,ii) = sign(x(ii,4:7)').*abs(q_dav(:,ii));
    error_dav(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_dav(:,ii)));
    q_svd(:,ii) = sign(x(ii,4:7)').*abs(q_svd(:,ii));
    error_svd(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_svd(:,ii)));
end

% Compute Mean Error and Plot Results
triad_err = mean(error_triad)*180/pi
dav_err = mean(error_dav)*180/pi
svd_err = mean(error_svd)*180/pi
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_triad(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Triad Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_triad(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_triad(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_triad(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Triad')
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_dav(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. Davenport Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_dav(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_dav(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_dav(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','Davenport')
figure
subplot(4,1,1),plot(t,q_true(:,1),t,q_svd(1,:),'r-'),ylabel('q_1'),xlabel('t'),title('Truth vs. SVD Quaternion')
subplot(4,1,2),plot(t,q_true(:,2),t,q_svd(2,:),'r-'),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3),t,q_svd(3,:),'r-'),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4),t,q_svd(4,:),'r-'),ylabel('q_4'),xlabel('t')
legend('Truth','SVD')


%% Static Estimation Monte Carlo Simulation

MC_len = 250;

for jj = 1:MC_len
    % Noisy Magnetometer Body-Frame Measurements
    rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
    rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
    rN = [rN1,rN2];
    [rB_noisy,~] = generateNoise(x,V_mag,1,rN);
    for ii = 1:length(x)
        rB = [rB_noisy(1:3,ii) rB_noisy(4:6,ii)];
        [q_triad(:,ii),~] = triad(rN,rB);
        [q_dav(:,ii),~] = qdaven(rN,rB,error_mag);
        [q_svd(:,ii),~] = qsvd(rN,rB,error_mag);
        q_triad(:,ii) = sign(x(ii,4:7)').*abs(q_triad(:,ii));
        error_triad(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_triad(:,ii)));
        q_dav(:,ii) = sign(x(ii,4:7)').*abs(q_dav(:,ii));
        error_dav(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_dav(:,ii)));
        q_svd(:,ii) = sign(x(ii,4:7)').*abs(q_svd(:,ii));
        error_svd(:,ii) = norm(quat2phi(qmult(qconj(x(ii,4:7)'))*q_svd(:,ii)));
    end
    triad_mean(jj) = mean(error_triad)*180/pi;
    dav_mean(jj) = mean(error_dav)*180/pi;
    svd_mean(jj) = mean(error_svd)*180/pi;
    
    if mod(jj,25) == 0
        fprintf('Completed %0.0f/%0.0f Simulations\n',jj,MC_len);
    end
    
end

% Plot Mean Error Results
figure
subplot(3,1,1),plot(1:1:jj,triad_mean,'.'),ylabel('Error (deg)')
title('TRIAD Method Error')
subplot(3,1,2),plot(1:1:jj,dav_mean,'.'),ylabel('Error (deg)')
title('Davenport q-method Method Error')
subplot(3,1,3),plot(1:1:jj,svd_mean,'.'),xlabel('Samples'),ylabel('Error (deg)')
title('SVD Method Error')

%% Recursive Attitude Estimation (MEKF)
clc,clear,close all
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')

% From Superspin
rho = [0;0;20.8612];

% Initial Attitude Conditions
q0 = [0;0;0;1]; % Identity Quaternion
om0 = [0.25;0.25;2.5]; % rad/s
x0 = [om0;q0];

% ODE45 Solver for Attitude Dynamics
dt = 0.01;
tspan = 0:dt:6;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)SkySat_Dynamics_Sim(t,x,mu,J,rho,0), tspan, x0, opts);

% Plot Quaternion
om_true = x(:,1:3); q_true = x(:,4:7);

% Noisy Magnetometer Body-Frame Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = generateNoise(x,V_mag,1,rN);

% Noisy Star Tracker Quaternion Measurements
[q_noisy,~] = generateNoise(x,V_st,2,[]);

% Noisy Gyro State Measurements (with Bias)
[om_noisy,bias] = generateNoise(x,W_gyro,3,[]);

% MEKF Noisy Measurement Histories
whist = om_noisy;
yhist = [q_noisy; rB_noisy];

% MEKF Initial Values
q0 = q_true(1,:)';
x0 = [q0;0;0;0]; % [quaternion; gyro bias]
P0 = (10*pi/180)^2*eye(6);
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

P0 = (pi)^2*eye(6);
W_MEKF = W_gyro*100000;
V_MEKF = blkdiag(V_st*10, V_mag*10, V_mag*10);

tic

% Run MEKF
[xhist, Phist] = mekf(x0, P0, W_MEKF, V_MEKF, rN, whist, yhist, dt);

% Calculate Error Quaternions
e = zeros(3,size(q_true,1));
for k = 1:size(q_true,1)
    e(:,k) = quat2phi(qmult(qconj(q_true(k,:)'))*xhist(1:4,k));
end

MEKF_time = toc

% Plot Attitude
figure,
subplot(4,1,1),hold on;...
    plot(t,xhist(1,:),'r--','LineWidth',2),plot(t,q_true(:,1),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_1');...
    title('MEKF Quaternion Estimation');legend('True', 'Estimated')
subplot(4,1,2),hold on;...
    plot(t,xhist(2,:),'r--','LineWidth',2),plot(t,q_true(:,2),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_2');
subplot(4,1,3),hold on;...
    plot(t,xhist(3,:),'r--','LineWidth',2),plot(t,q_true(:,3),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_3');
subplot(4,1,4),hold on;...
    plot(t,xhist(4,:),'r--','LineWidth',2),plot(t,q_true(:,4),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_4');

% Plot Errors
figure;
subplot(3,1,1);
plot(t,(180/pi)*e(1,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');ylim([-1,1]);
subplot(3,1,2);
plot(t,(180/pi)*e(2,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('Degrees');ylim([-1,1]);
subplot(3,1,3);
plot(t,(180/pi)*e(3,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
xlabel('Time (s)');ylim([-1,1]);

figure;
subplot(3,1,1);
plot(t,xhist(5,:)-bias(1,:)); hold on
plot(t,2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error'); 
subplot(3,1,2);
plot(t,xhist(6,:)-bias(2,:)); hold on
plot(t,2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(5,5,:))),'r');
subplot(3,1,3);
plot(t,xhist(7,:)-bias(3,:)); hold on
plot(t,2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(6,6,:))),'r');
xlabel('Time (s)');

%% Recursive Attitude Estimation (UKF)

clc,clear
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')

tic 

% From Superspin
rho = [0;0;20.8612];

% Initial Attitude Conditions
q0 = [0;0;0;1];
om0 = [0.25;0.25;2.5];
x0 = [om0;q0];

% ODE45 Solver for Attitude Dynamics
dt = 0.01;
tspan = 0:dt:6;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)SkySat_Dynamics_Sim(t,x,mu,J,rho,0), tspan, x0, opts);

% Plot Quaternion
om_true = x(:,1:3); q_true = x(:,4:7);

% Noisy Magnetometer Body-Frame Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = generateNoise(x,V_mag,1,rN);

% Noisy Star Tracker Quaternion Measurements
[q_noisy,~] = generateNoise(x,V_st,2,[]);

% Noisy Gyro State Measurements (with Bias)
[om_noisy,bias] = generateNoise(x,W_gyro,3,[]);

% MEKF Noisy Measurement Histories
whist = om_noisy;
yhist = [q_noisy; rB_noisy];

% UKF Initial Values
q0 = q_true(1,:)';
x0 = [q0;0;0;0]; % [quaternion; gyro bias]
P0 = (pi)^2*eye(6);
W_MEKF = W_gyro*10000;
V_MEKF = blkdiag(V_st*10, V_mag*10, V_mag*10);

tic

% Run UKF
[xhist, Phist] = ukf(x0, P0, W_MEKF, V_MEKF, rN, whist, yhist, dt);

% Calculate Error Quaternions
e = zeros(3,size(q_true,1));
for k = 1:size(q_true,1)
    e(:,k) = quat2phi(qmult(qconj(q_true(k,:)'))*xhist(1:4,k));
end

UKF_time = toc

% Plot Attitude
figure,
subplot(4,1,1),hold on;...
    plot(t,xhist(1,:),'r--','LineWidth',2),plot(t,q_true(:,1),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_1');...
    title('UKF Quaternion Estimation');legend('True', 'Estimated')
subplot(4,1,2),hold on;...
    plot(t,xhist(2,:),'r--','LineWidth',2),plot(t,q_true(:,2),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_2');
subplot(4,1,3),hold on;...
    plot(t,xhist(3,:),'r--','LineWidth',2),plot(t,q_true(:,3),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_3');
subplot(4,1,4),hold on;...
    plot(t,xhist(4,:),'r--','LineWidth',2),plot(t,q_true(:,4),'LineWidth',1);...
    hold off, xlabel('Time (s)'),ylabel('q_4');

% Plot Errors
figure;
subplot(3,1,1);
plot(t,(180/pi)*e(1,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(1,1,:)))/40,'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(1,1,:)))/40,'r');
title('Attitude Error');ylim([-1,1]);
subplot(3,1,2);
plot(t,(180/pi)*e(2,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(2,2,:)))/40,'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(2,2,:)))/40,'r');
ylabel('Degrees');ylim([-1,1]);
subplot(3,1,3);
plot(t,(180/pi)*e(3,:)); hold on
plot(t,2*(180/pi)*sqrt(squeeze(Phist(3,3,:)))/40,'r');
plot(t,-2*(180/pi)*sqrt(squeeze(Phist(3,3,:)))/40,'r');
xlabel('Time (s)');ylim([-1,1]);

figure;
subplot(3,1,1);
plot(t,xhist(5,:)-bias(1,:)); hold on
plot(t,2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error');ylim([-.5,.5]);
subplot(3,1,2);
plot(t,xhist(6,:)-bias(2,:)); hold on
plot(t,2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(5,5,:))),'r');ylim([-.5,.5]);
subplot(3,1,3);
plot(t,xhist(7,:)-bias(3,:)); hold on
plot(t,2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(t,-2*sqrt(squeeze(Phist(6,6,:))),'r');
xlabel('Time (s)');ylim([-.5,.5]);

%% SkySat1 Dynamics with Environmental Disturbances
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\geometry')

% From Superspin
rho = [0;0;20.8612];

% Initial Conditions
[r_eci, v_eci] = OE2ECI(a, e, deg2rad(inc), deg2rad(RAAN), deg2rad(w), anom, mu);
q0 = [0;0;0;1]; % Identity Quaternion
om0 = [0;0;0]; % rad/s, Start From Rest
x0 = [om0;q0;r_eci;v_eci];

% ODE45 Solver for Full SkySat1 Dynamics
dt = 10;
tspan = 0:dt:15*T;
tol = 1e-6; opts  = odeset('reltol', tol, 'abstol', tol);
[t, x] = ode45(@(t,x)SkySat_Dynamics_Sim(t,x,mu,J,rho,params,1), tspan, x0, opts);

% % Plot Orbit
% figure, hold on
% plot3(x(:,8),x(:,9),x(:,10),'r','LineWidth',2),axis equal
% xlabel('[km]'),ylabel('[km]'),zlabel('[km]')
% title('SkySat1 Orbit Dynamics')
% earthPlot(1), axis equal
% hold off

% % Plot Angular Velocity
% om_true = x(:,1:3);
% figure
% subplot(3,1,1),plot(t,om_true(:,1)),ylabel('\omega_x'),xlabel('t'),title('SkySat1 Angular Velocity with Disturbances')
% subplot(3,1,2),plot(t,om_true(:,2)),ylabel('\omega_y'),xlabel('t')
% subplot(3,1,3),plot(t,om_true(:,3)),ylabel('\omega_z'),xlabel('t')

% Plot Quaternion
q_true = x(:,4:7);
figure
subplot(4,1,1),plot(t,q_true(:,1)),ylabel('q_1'),xlabel('t'),title('SkySat1 Quaternion Dynamics with Disturbances')
subplot(4,1,2),plot(t,q_true(:,2)),ylabel('q_2'),xlabel('t')
subplot(4,1,3),plot(t,q_true(:,3)),ylabel('q_3'),xlabel('t')
subplot(4,1,4),plot(t,q_true(:,4)),ylabel('q_4'),xlabel('t')

% Plot Orbital Decay
orb_alt = vecnorm(x(:,8:10)');
FO = fit(t, orb_alt', 'poly1')
bestfit = -1.259e-05*t + 6950;

figure;
hold on; plot(t/T,orb_alt),plot(t/T,bestfit,'r','LineWidth',2),...
    title('Orbital Decay Over 15 Periods'),xlabel('Orbit Count'),...
    ylabel('Radial Distance from Earth (km)');grid on;

total_decay = bestfit(1)-bestfit(end)

% Plot Angular Momentum
for ii = 1:size(x,1)
    dt(ii) = norm(J*x(ii,1:3)');
end
figure;
plot(t./T,dt,'b','LineWidth',1)
title('Angular Momentum')
xlabel('Orbit Count')
ylabel('Angular Momentum (kg-m\^2/sec)')


%% Recover Air Drag and Torque History
for ii = 1:size(x,1)
    Q = quat2Q(x(ii,4:7)');
    [Drag(:,ii),Torque(:,ii)] = SkySat_Drag_Sim(x(ii,8:10)',x(ii,11:13)',Q,S,n,c);
    Tau_g(:,ii) = cross(3*mu/(dot(x(ii,8:10)',x(ii,8:10)')^(5/2))*x(ii,8:10)', Q*J/(1000^2)*x(ii,8:10)'); % [N-km]
    Tau_g(:,ii) = Q'*Tau_g(:,ii)*1000; % Rotate into Body Frame and Convert Units, [N-m]
end

max_drag = max(vecnorm(Drag));
max_dragtorque = max(vecnorm(Torque))
max_gravtorque = max(vecnorm(Tau_g))
average_torque = mean(vecnorm(Torque + Tau_g));

%% Simulate LQR-Controlled Attitude Dynamics with Disturbances and MEKF
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')
load('Parameters\geometry')

% Initial Conditions
[r_eci, v_eci] = OE2ECI(a, e, deg2rad(inc), deg2rad(RAAN), deg2rad(w), anom, mu);
qstar = [1/sqrt(2); 1/sqrt(2); 0; 0]; % Desired Quaternion
phi_rand = rand(3,1);
phi0 = (pi/2)*phi_rand/norm(phi_rand); % Normalize and Set 90° Initial Error
%phi0 = (pi/4)*phi_rand/norm(phi_rand); % Normalize and Set 45° Initial Error
%phi0 = (pi/8)*phi_rand/norm(phi_rand); % Normalize and Set 22.5° Initial Error
%phi0 = (pi/16)*phi_rand/norm(phi_rand); % Normalize and Set 11.25° Initial Error
q_err = phi2quat(phi0); q(:,1) = qmult(q_err)*qstar;
om0 = [0; 0; 0]; % rad/s, Zero Initial Angular Velocity
x0 = [phi0; om0; r_eci; v_eci];

% Compute LQR Gains about Linearized Dynamics
A = [zeros(3) eye(3); zeros(3) zeros(3)];
B = [zeros(3); -inv(J)];
Q = (1/(10000)^2)*eye(6); R = 1*eye(size(B,2)); % Q is Tuned Based on Assumed Initial Error
[K,~,~] = lqr(A,B,Q,R);

% Initial Noisy Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = generateNoise([0,0,0,q(:,1)'],V_mag,1,rN); % Magnetometer
[q_noisy,~] = generateNoise([0,0,0,q(:,1)'],V_st,2,[]); % Star Tracker
[om_noisy,bias] = generateNoise([0,0,0,q(:,1)'],W_gyro,3,[]); % Gyro
    
% MEKF Initial Values
xhist(:,1) = [q(:,1);0;0;0]; % [quaternion; gyro bias]
Phist(:,:,1) = (10*pi/180)^2*eye(6);
whist(:,1) = om_noisy;
yhist(:,1) = [q_noisy; rB_noisy];
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

% Runge-Kutta Integrator for Controlled Attitude Dynamics
dt = 1; tspan = (0:dt:4*T)/60; % Minutes
x = zeros(length(x0),length(tspan)); x(:,1) = x0; 
theta = zeros(1,length(tspan)); theta(1) = norm(phi0)*180/pi;
u = zeros(3,length(tspan)); xhat = x0(1:6);
for ii = 2:length(tspan)
    % Feedback Control Law
    u(:,ii) = -K*xhat; % u = rho_dot
    % RK4 Step
    k1 = dt*ODEstep([],x(:,ii-1),u(:,ii),A,B,q(:,ii-1),params);
    k2 = dt*ODEstep([],x(:,ii-1)+k1/2,u(:,ii),A,B,q(:,ii-1),params);
    k3 = dt*ODEstep([],x(:,ii-1)+k2/2,u(:,ii),A,B,q(:,ii-1),params);
    k4 = dt*ODEstep([],x(:,ii-1)+k3,u(:,ii),A,B,q(:,ii-1),params);
    x(:,ii) = x(:,ii-1) + (k1+2*k2+2*k3+k4)/6;
    % Compute Angle Error, Quaternion, and Actuator Speeds
    theta(ii) = norm(x(1:3,ii))*180/pi;
    q(:,ii) = qmult(phi2quat(x(1:3,ii)))*qstar;
    inputs(:,ii) = pinv(Bw)*u(:,ii);
    % Noisy Measurements at a Single Time Step
    [rB_noisy,~] = generateNoise([0,0,0,q(:,ii)'],V_mag,1,rN); % Magnetometer
    [q_noisy,~] = generateNoise([0,0,0,q(:,ii)'],V_st,2,[]); % Star Tracker
    [om_noisy,bias] = generateNoise([0,0,0,q(:,ii)'],W_gyro,3,[]); % Gyro

    % MEKF Noisy Measurement Histories
    whist(:,ii) = om_noisy;
    yhist(:,ii) = [q_noisy; rB_noisy];

    % Estimate Attitude with MEKF at a Single Time Step
    [xhist(:,ii), Phist(:,:,ii)] = mekf_step(xhist(:,ii-1), Phist(:,:,ii-1), W_MEKF, V_MEKF, rN, whist(:,ii), yhist(:,ii), dt);
    xhat = [quat2phi(qmult(qconj(qstar))*xhist(1:4,ii)); x(4:6,ii)];
    xhat(1:3,:) = [xhat(2,:); xhat(1,:); -xhat(3,:)];
end
% Check Torque and Speed Values; Rxn Wheel Speed Limit: 3200 RPM(2pi/60) = 335 rad/s
Max_Torque1 = max(inputs(1,:)) % [N-m] Torque Limit per Wheel: 0.82 N-m
Max_Torque2 = max(inputs(2,:))
Max_Torque3 = max(inputs(3,:))
Max_Torque4 = max(inputs(4,:))

% Plot Attitude Error
figure, plot(tspan,theta),ylabel('Degrees'),xlabel('Time (min)')
title('Attitude Error')

% Plot Angular Velocity
om = x(4:6,:); figure
subplot(3,1,1),plot(tspan,om(1,:)),ylabel('\omega_x'),xlabel('Time (s)')
title('SkySat1 Angular Velocity with Disturbances')
subplot(3,1,2),plot(tspan,om(2,:)),ylabel('\omega_y'),xlabel('Time (s)')
subplot(3,1,3),plot(tspan,om(3,:)),ylabel('\omega_z'),xlabel('Time (s)')

% Plot Axis-Angle Error
phi = x(1:3,:); figure
subplot(3,1,1),plot(tspan,phi(1,:)),ylabel('\phi_1'),xlabel('Time (s)')
title('Axis-Angle Error')
subplot(3,1,2),plot(tspan,phi(2,:)),ylabel('\phi_2'),xlabel('Time (s)')
subplot(3,1,3),plot(tspan,phi(3,:)),ylabel('\phi_3'),xlabel('Time (s)')

% Plot True vs. Estimated Quaternion
figure,
subplot(4,1,1),plot(tspan,q(1,:),tspan,xhist(1,:),'r-'),xlabel('Time (s)'),ylabel('q_1')
title('SkySat1 Quaternion Dynamics'); legend('True', 'Estimated')
subplot(4,1,2),plot(tspan,q(2,:),tspan,xhist(2,:),'r-'),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(tspan,q(3,:),tspan,xhist(3,:),'r-'),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(tspan,q(4,:),tspan,xhist(4,:),'r-'),xlabel('Time (s)'),ylabel('q_4')

% Calculate and Plot Quaternion Estimate Error
e = zeros(3,size(q,2));
for k = 1:size(q,2)
    e(:,k) = quat2phi(qmult(qconj(q(:,k)))*xhist(1:4,k));
end
figure;
subplot(3,1,1);
plot(tspan,(180/pi)*e(1,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');
subplot(3,1,2);
plot(tspan,(180/pi)*e(2,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('Degrees');
subplot(3,1,3);
plot(tspan,(180/pi)*e(3,:)); axis([0 tspan(end) -0.25 0.25]); hold on
plot(tspan,2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(tspan,-2*(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
xlabel('Time (s)')

% Plot Reaction Wheel Speeds
figure
subplot(4,1,1),plot(tspan,inputs(1,:)),ylabel('\tau_{rxn}^1'),xlabel('t'),title('RXN Wheel Torque')
subplot(4,1,2),plot(tspan,inputs(2,:)),ylabel('\tau_{rxn}^2'),xlabel('t')
subplot(4,1,3),plot(tspan,inputs(3,:)),ylabel('\tau_{rxn}^3'),xlabel('t')
subplot(4,1,4),plot(tspan,inputs(4,:)),ylabel('\tau_{rxn}^4'),xlabel('t')

%% Eigen-Axis Slew
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\geometry')

% Nominal Versine Function
theta_f = 179*pi/180; % Radians
alpha = pi/1600; dt = 0.1; tspan = 0:dt:pi/alpha;
theta = 0.5*theta_f*(1 - cos(alpha*tspan)); % Versine Function
theta_dot = 0.5*alpha*theta_f*sin(alpha*tspan);
theta_ddot = 0.5*(alpha^2)*theta_f*cos(alpha*tspan);

% Nominal Reference Trajectory
phi_rand = rand(3,1);
r = phi_rand/norm(phi_rand); % Unit Vector, Axis of Rotation
phi_traj = theta.*r;
om_traj = theta_dot.*r;
om_dot_traj = theta_ddot.*r;
for ii = 1:length(phi_traj)
    q_traj(:,ii) = phi2quat(phi_traj(:,ii));
end

% Solve for Nominal Reaction Wheel Inputs with 4th Order Runge Kutta
rho_traj(:,1) = [0;0;0]; rho_dot(:,1) = [0;0;0];
rho_dots = @(t,rho,om,om_dot) -J*om_dot - cross(om,J*om + rho);
for ii = 2:length(tspan)
    % RK4 Step
    k1 = rho_dots([],rho_traj(:,ii-1),om_traj(:,ii-1),om_dot_traj(:,ii-1));
    k2 = rho_dots([],rho_traj(:,ii-1)+k1/2,om_traj(:,ii-1),om_dot_traj(:,ii-1));
    k3 = rho_dots([],rho_traj(:,ii-1)+k2/2,om_traj(:,ii-1),om_dot_traj(:,ii-1));
    k4 = rho_dots([],rho_traj(:,ii-1)+k3,om_traj(:,ii-1),om_dot_traj(:,ii-1));
    rho_dot(:,ii) = (k1+2*k2+2*k3+k4)/6; % Nominal Input Torque
    rho_traj(:,ii) = rho_traj(:,ii-1) + dt*rho_dot(:,ii);
    % Map into Inputs and Speed
    inputs_traj(:,ii) = pinv(Bw)*rho_dot(:,ii); % Inverse Dynamics for Wheel Torque
    speeds_traj(:,ii) = pinv(Bw)*rho_traj(:,ii); % Inverse Dynamics for Wheel Speed
end

% Linearize and Discretize About Nominal Trajectory for LQR Gains
Q1 = ((1/(pi))^2)*ones(1,3); Q2 = ((1/(pi))^2)*ones(1,3); Q3 = ((1/(pi))^2)*ones(1,3);
Q = diag([Q1 Q2 Q3]); S = 1000*Q;
R = 1*eye(3);
for ii = 2:length(tspan)
    Ac = [-hat(om_traj(:,ii-1)) eye(3) zeros(3);
          zeros(3) -inv(J)*(hat(om_traj(:,ii-1))*J-hat(J*om_traj(:,ii-1)+rho_traj(:,ii-1))) -inv(J)*hat(om_traj(:,ii-1));
          zeros(3) zeros(3) zeros(3)];
    Bc = [zeros(3); -inv(J); eye(3)]; % u = rho_dot
    At = eye(9) + Ac*dt; Bt = Bc*dt; % Short Time-Step Approximation
    [K(:,:,ii),S] = tvlqr(At,Bt,Q,R,S); % Compute Time-History of LQR Gains
end

save('nom_traj','phi_traj','om_traj','q_traj','rho_traj','rho_dot','inputs_traj','speeds_traj','K')
disp('Saved Parameters.')

% Time-Varying LQR for Attitude Tracking with MEKF
clc,clearvars
load('Parameters\orbital')
load('Parameters\mass')
load('Parameters\sensor')
load('Parameters\geometry')
load('nom_traj') % Nominal Values: phi_traj, om_traj, q_traj, rho, rho_dot

% Initial State and Perturbation Conditions
alpha = pi/1600; dt = 0.1; tspan = 0:dt:pi/alpha;
[r_eci, v_eci] = OE2ECI(a, e, deg2rad(inc), deg2rad(RAAN), deg2rad(w), anom, mu);
dphi(:,1) = [0; 0; 5*pi/180]; % Normalize and Set Initial Axis-Angle Perturbation
dom(:,1) = [0;0;0]; % rad/s, Zero Initial Angular Velocity Perturbation
drho(:,1) = [0;0;0]; % Zero Initial Wheel Speed Perturbation
dq(:,1) = phi2quat(dphi(:,1)); % Initial Quaternion Perturbation
dx(:,1) = [dphi(:,1); dom(:,1); drho(:,1)]; 

phi(:,1) = phi_traj(:,1) + dphi(:,1);
om(:,1) = om_traj(:,1) + dom(:,1);
rho(:,1) = rho_traj(:,1) + drho(:,1);
q(:,1) = qmult(q_traj(:,1))*dq(:,1);
x(:,1) = [q(:,1); om(:,1); rho(:,1); r_eci; v_eci]; % [q; om; rho; r; v]

% Initial Noisy Measurements
rN1 = 2*rand(3,1)-1; rN1 = rN1/norm(rN1);
rN2 = 2*rand(3,1)-1; rN2 = cross(rN1,rN2); rN2 = rN2/norm(rN2);
rN = [rN1,rN2];
[rB_noisy,~] = generateNoise([0,0,0,q(:,1)'],V_mag,1,rN); % Magnetometer
[q_noisy,~] = generateNoise([0,0,0,q(:,1)'],V_st,2,[]); % Star Tracker
[om_noisy,bias] = generateNoise([0,0,0,q(:,1)'],W_gyro,3,[]); % Gyro
    
% MEKF Initial Values
xhist(:,1) = [q(:,1);0;0;0]; % [quaternion; gyro bias]
Phist(:,:,1) = (10*pi/180)^2*eye(6);
whist(:,1) = om_noisy;
yhist(:,1) = [q_noisy; rB_noisy];
W_MEKF = W_gyro;
V_MEKF = blkdiag(V_st, V_mag, V_mag);

% Runge-Kutta Integrator for Controlled Attitude Tracking
u(:,1) = rho_dot(:,1);
for ii = 2:length(tspan)        
    % Compute Gains and Feedback Control Law
    du = -K(:,:,ii)*dx;
    u(:,ii) = rho_dot(:,ii) + du; % rho_dot = nominal control 
    
    % Discrete RK4 Step with Discrete Error Dynamics and Continuous Orbit Dynamics
    k1 = ODEtrack([],x(:,ii-1),u(:,ii-1),Bw,params);
    k2 = ODEtrack([],x(:,ii-1)+k1/2,u(:,ii-1),Bw,params);
    k3 = ODEtrack([],x(:,ii-1)+k2/2,u(:,ii-1),Bw,params);
    k4 = ODEtrack([],x(:,ii-1)+k3,u(:,ii-1),Bw,params);
    x(:,ii) = x(:,ii-1) + dt*(k1+2*k2+2*k3+k4)/6; % Continuous Orbit Dynamics
    
    % Parse Propagated State [q; om; rho]
    q(:,ii) = x(1:4,ii)/norm(x(1:4,ii));
    phi(:,ii) = quat2phi(q(:,ii));  
    om(:,ii) = x(5:7,ii);
    rho(:,ii) = x(8:10,ii);
    
    % Compute Angle Error and Actuator Speeds
    theta(ii) = norm(phi(:,ii))*180/pi;
    inputs(:,ii) = pinv(Bw)*u(:,ii); % Scalar Torques
    speeds(:,ii) = pinv(Bw)*x(8:10,ii);
    
    % Noisy Measurements at a Single Time Step
    [rB_noisy,~] = generateNoise([0,0,0,q(:,ii)'],V_mag,1,rN); % Magnetometer
    [q_noisy,~] = generateNoise([0,0,0,q(:,ii)'],V_st,2,[]); % Star Tracker
    [om_noisy,bias] = generateNoise([0,0,0,q(:,ii)'],W_gyro,3,[]); % Gyro
    whist(:,ii) = om_noisy;
    yhist(:,ii) = [q_noisy; rB_noisy];
    
    % Estimate Attitude with MEKF at a Single Time Step
    [xhist(:,ii),Phist(:,:,ii)] = mekf_step(xhist(:,ii-1),Phist(:,:,ii-1),W_MEKF,V_MEKF,rN,whist(:,ii),yhist(:,ii),dt);
    dx = [quat2phi(qmult(qconj(q_traj(:,ii)))*xhist(1:4,ii)); x(5:7,ii) - om_traj(:,ii); x(8:10,ii) - rho_traj(:,ii)];
    %dx(1:3,:) = [-dx(2,:); dx(1,:); -dx(3,:)]; % Sign Correction
end

% Compute Max Torques and Speeds
Max_Torque1 = max(abs(inputs(1,:))) % [N-m] Torque Limit for Each Rxn Wheels: 0.82 N-m
Max_Torque2 = max(abs(inputs(2,:)))
Max_Torque3 = max(abs(inputs(3,:)))
Max_Torque4 = max(abs(inputs(4,:)))

Max_Speed1 = max(abs(speeds(1,:))) % [rad/s] Rxn Wheel Speed Limit: 3200 RPM(2pi/60) = 335 rad/s
Max_Speed2 = max(abs(speeds(2,:)))
Max_Speed3 = max(abs(speeds(3,:)))
Max_Speed4 = max(abs(speeds(4,:)))

% Plot Nominal vs. Closed-Loop Results
figure
subplot(3,1,1),plot(tspan,om_traj(1,:),tspan,om(1,:)),xlabel('Time (s)'),ylabel('\omega_x')
title('Nominal vs. Closed-Loop Angular Velocity')
subplot(3,1,2),plot(tspan,om_traj(2,:),tspan,om(2,:)),xlabel('Time (s)'),ylabel('\omega_y')
subplot(3,1,3),plot(tspan,om_traj(3,:),tspan,om(3,:)),xlabel('Time (s)'),ylabel('\omega_z')
legend('Nominal','Closed-Loop','Location','Best')

figure
subplot(4,1,1),plot(tspan,q_traj(1,:),tspan,q(1,:)),xlabel('Time (s)'),ylabel('q_1')
title('Nominal vs. Closed-Loop Quaternion')
subplot(4,1,2),plot(tspan,q_traj(2,:),tspan,q(2,:)),xlabel('Time (s)'),ylabel('q_2')
subplot(4,1,3),plot(tspan,q_traj(3,:),tspan,q(3,:)),xlabel('Time (s)'),ylabel('q_3')
subplot(4,1,4),plot(tspan,q_traj(4,:),tspan,q(4,:)),xlabel('Time (s)'),ylabel('q_4')
legend('Nominal','Closed-Loop','Location','Best')

figure
plot(tspan,vecnorm(phi_traj)*180/pi,tspan,theta),xlabel('Time (s)'),ylabel('\theta')
title('Nominal vs. Closed-Loop Angle Error')
legend('Nominal','Closed-Loop','Location','Best')

figure
subplot(3,1,1),plot(tspan,phi_traj(1,:),tspan,phi(1,:)),xlabel('Time (s)'),ylabel('\phi_1')
title('Nominal vs. Closed-Loop Axis-Angle')
subplot(3,1,2),plot(tspan,phi_traj(2,:),tspan,phi(2,:)),xlabel('Time (s)'),ylabel('\phi_2')
subplot(3,1,3),plot(tspan,phi_traj(3,:),tspan,phi(3,:)),xlabel('Time (s)'),ylabel('\phi_3')
legend('Nominal','Closed-Loop','Location','Best')

figure
subplot(4,1,1),plot(tspan,inputs_traj(1,:),tspan,inputs(1,:)),ylabel('\tau_{rxn}^1'),xlabel('t')
title('Reaction Wheel Torques') % [N-m]
subplot(4,1,2),plot(tspan,inputs_traj(2,:),tspan,inputs(2,:)),ylabel('\tau_{rxn}^2'),xlabel('t')
subplot(4,1,3),plot(tspan,inputs_traj(3,:),tspan,inputs(3,:)),ylabel('\tau_{rxn}^3'),xlabel('t')
subplot(4,1,4),plot(tspan,inputs_traj(4,:),tspan,inputs(4,:)),ylabel('\tau_{rxn}^4'),xlabel('t')
legend('Nominal','Closed-Loop','Location','Best')

figure
subplot(4,1,1),plot(tspan,speeds_traj(1,:),tspan,speeds(1,:)),ylabel('\omega_{rxn}^1'),xlabel('t')
title('Reaction Wheel Speeds') % [rad/s]
subplot(4,1,2),plot(tspan,speeds_traj(1,:),tspan,speeds(2,:)),ylabel('\omega_{rxn}^2'),xlabel('t')
subplot(4,1,3),plot(tspan,speeds_traj(1,:),tspan,speeds(3,:)),ylabel('\omega_{rxn}^3'),xlabel('t')
subplot(4,1,4),plot(tspan,speeds_traj(1,:),tspan,speeds(4,:)),ylabel('\omega_{rxn}^4'),xlabel('t')
legend('Nominal','Closed-Loop','Location','Best')
    