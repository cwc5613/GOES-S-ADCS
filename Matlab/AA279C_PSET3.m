clear;
clc;
close all;

[a, e, inc, RAAN, w, M0, E, anom, revs_per_day] =...
    TLE_Reader('SKYSAT1_TLE.txt');

load 'AA279C_ps2_data.mat';

b_R_P = [Ix;Iy;Iz];

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
tol = 1e-6;

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci; v_eci];

stopTime = 0.1*60;
stepSize = 1/2000;
tspan = [0:stepSize:stopTime];
opts  = odeset('reltol', tol, 'abstol', tol);

[t,y] = ode45(@(t,y) orb_prop(y,mu), tspan, initCond, opts);

% figure
% hold on
% plot3(y(:,1),y(:,2),y(:,3),'c','LineWidth',1)
% earthPlot(1)
% axis equal
% hold off
% xlabel('[km]')
% ylabel('[km]')
% zlabel('[km]')

% [t2,y2] = ode45(@(t,y) J2_prop(y,mu,J2,Re), tspan, initCond, opts);

% figure
% hold on
% plot3(y2(:,1),y2(:,2),y2(:,3),'c','LineWidth',1)
% earthPlot(1)
% axis equal
% hold off
% xlabel('[km]')
% ylabel('[km]')
% zlabel('[km]')

%% Quaternions

% Initial Attitude Conditions
Q0 = [0 -1 0;1 0 0;0 0 1];
q0 = Q2quat(Q0);
om0 = [0;1;10.7];
x0 = [r_eci;v_eci;om0;q0];

% ODE45 Solver for Orbital Dynamics
[t, y3] = ode45(@(t,y3)HSTdynamics(mu,J,rho,y3), tspan, x0, opts);
qtrue = y3(:,10:13);

%Star Tracker Measurements
ST_error = [18;18;102].*0.00000484813; %radians

% Star Tracker Error
for i = 1:size(y3,1)
    ST_phi = (randn(1,3)*chol(diag(ST_error)))';
    q_error = [ST_phi/2;1-(1/8)*(ST_phi'*ST_phi)];
    q_error = q_error/norm(q_error);
    q_noise(i,:) = (qmult(qtrue(i,:)')*q_error)';
end

%Gyroscope Measurements
GYRO_RND = 0.018*pi/180; %deg/s/sqrt(Hz)
GYRO_ARW = 0.002*pi/180*60; %deg/sqrt(Hr)
GYRO_error = GYRO_ARW^2*eye(3);

%GPS Measurements
GPS_var = 0.02*pi/180*60; %deg/sqrt(Hr)
GPS_error = GPS_var^2*eye(3);

rn1 = 2*rand(3,1)-1;
rn1 = rn1/norm(rn1);
rn2 = 2*rand(3,1)-1;
rn2 = rn2/norm(rn2);
rn3 = cross(rn1,rn2);
rn3 = rn3/norm(rn3);
rn2 = cross(rn1,rn3);
rn2 = rn2/norm(rn2);

for i = 1:size(qtrue,1)
    GYRO_err = (randn(1,3)*chol(01*GYRO_error))';
    GPS_err = (randn(1,3)*chol(01*GPS_error))';
    Q = q2Q(qtrue(i,:))';
    rb1(:,i) = Q*rn1 + GYRO_err;
    rb1(:,i) = rb1(:,i)/norm(rb1(:,i));
    rb3(:,i) = Q*rn3 + GYRO_err;
    rb3(:,i) = rb3(:,i)/norm(rb3(:,i));
    rb2(:,i) = Q*rn2 + GPS_err;
    rb2(:,i) = rb2(:,i)/norm(rb2(:,i));
    [q1(:,i),Q1] = svd_method(rn1,rn2,rn3,rb1(:,i),rb2(:,i),rb3(:,i),GYRO_ARW,GPS_var,GYRO_ARW);
    [q2(:,i),Q2] = davenport_method(rn1,rn2,rn3,rb1(:,i),rb2(:,i),rb3(:,i),GYRO_ARW,GPS_var,GYRO_ARW);
    [q3(:,i),Q3] = triad(rn1,rn2,rb1(:,i),rb2(:,i));
    
    q1(:,i) = sign(qtrue(i,:))'.*abs(q1(:,i));
    q2(:,i) = sign(qtrue(i,:))'.*abs(q2(:,i));
    q3(:,i) = sign(qtrue(i,:))'.*abs(q3(:,i));
    
end

% figure
% subplot(3,1,1),plot(t,rb1(1,:),'LineWidth',2),ylabel('q_1'),xlabel('Time (sec)'),title('Perturbed Rotor-Stabilized Quaternion Components vs Time')
% subplot(3,1,2),plot(t,rb1(2,:),'LineWidth',2),ylabel('q_2'),xlabel('Time (sec)')
% subplot(3,1,3),plot(t,rb1(3,:),'LineWidth',2),ylabel('q_3'),xlabel('Time (sec)')

% figure
% subplot(3,1,1),plot(t,rb2(1,:),'LineWidth',2),ylabel('q_1'),xlabel('Time (sec)'),title('Perturbed Rotor-Stabilized Quaternion Components vs Time')
% subplot(3,1,2),plot(t,rb2(2,:),'LineWidth',2),ylabel('q_2'),xlabel('Time (sec)')
% subplot(3,1,3),plot(t,rb2(3,:),'LineWidth',2),ylabel('q_3'),xlabel('Time (sec)')

% % Plot Orbit
% figure, hold on
% plot3(y3(:,1),y3(:,2),y3(:,3),'c','LineWidth',2)
% xlabel('[km]')
% ylabel('[km]')
% zlabel('[km]')
% title('Perturbed Rotor-Stabilized Orbit')
% earthPlot(1)
% axis equal
% hold off

% Plot Quaternion
figure
subplot(4,1,1),plot(t,qtrue(:,1),'LineWidth',2),ylabel('q_1'),xlabel('Time (sec)'),title('Perturbed Rotor-Stabilized Quaternion Components vs Time')
subplot(4,1,2),plot(t,qtrue(:,2),'LineWidth',2),ylabel('q_2'),xlabel('Time (sec)')
subplot(4,1,3),plot(t,qtrue(:,3),'LineWidth',2),ylabel('q_3'),xlabel('Time (sec)')
subplot(4,1,4),plot(t,qtrue(:,4),'LineWidth',2),ylabel('q_4'),xlabel('Time (sec)')

% Plot Quaternion with Noise
figure
subplot(4,1,1),plot(t,q_noise(:,1),'LineWidth',2),ylabel('q_1'),xlabel('Time (sec)'),title('Perturbed Rotor-Stabilized Quaternion Components vs Time')
subplot(4,1,2),plot(t,q_noise(:,2),'LineWidth',2),ylabel('q_2'),xlabel('Time (sec)')
subplot(4,1,3),plot(t,q_noise(:,3),'LineWidth',2),ylabel('q_3'),xlabel('Time (sec)')
subplot(4,1,4),plot(t,q_noise(:,4),'LineWidth',2),ylabel('q_4'),xlabel('Time (sec)')

figure
subplot(4,1,1),hold on,plot(t,qtrue(:,1),'LineWidth',2),plot(t,q1(1,:),'r','LineWidth',1),ylabel('q_1'),xlabel('Time (sec)'),title('SVD Quaternion Components vs Time'),hold off
legend('true','svd')
subplot(4,1,2),hold on,plot(t,qtrue(:,2),'LineWidth',2),plot(t,q1(2,:),'r','LineWidth',1),ylabel('q_2'),xlabel('Time (sec)'),hold off
subplot(4,1,3),hold on,plot(t,qtrue(:,3),'LineWidth',2),plot(t,q1(3,:),'r','LineWidth',1),ylabel('q_3'),xlabel('Time (sec)'),hold off
subplot(4,1,4),hold on,plot(t,qtrue(:,4),'LineWidth',2),plot(t,q1(4,:),'r','LineWidth',1),ylabel('q_4'),xlabel('Time (sec)'),hold off

figure
subplot(4,1,1),hold on,plot(t,qtrue(:,1),'LineWidth',2),plot(t,q2(1,:),'r','LineWidth',1),ylabel('q_1'),xlabel('Time (sec)'),title('Davenport Quaternion Components vs Time')
legend('true','dav')
subplot(4,1,2),hold on,plot(t,qtrue(:,2),'LineWidth',2),plot(t,q2(2,:),'r','LineWidth',1),ylabel('q_2'),xlabel('Time (sec)')
subplot(4,1,3),hold on,plot(t,qtrue(:,3),'LineWidth',2),plot(t,q2(3,:),'r','LineWidth',1),ylabel('q_3'),xlabel('Time (sec)')
subplot(4,1,4),hold on,plot(t,qtrue(:,4),'LineWidth',2),plot(t,q2(4,:),'r','LineWidth',1),ylabel('q_4'),xlabel('Time (sec)')

figure
subplot(4,1,1),hold on,plot(t,qtrue(:,1),'LineWidth',2),plot(t,q3(1,:),'r','LineWidth',1),ylabel('q_1'),xlabel('Time (sec)'),title('TRIAD Quaternion Components vs Time')
legend('true','triad')
subplot(4,1,2),hold on,plot(t,qtrue(:,2),'LineWidth',2),plot(t,q3(2,:),'r','LineWidth',1),ylabel('q_2'),xlabel('Time (sec)')
subplot(4,1,3),hold on,plot(t,qtrue(:,3),'LineWidth',2),plot(t,q3(3,:),'r','LineWidth',1),ylabel('q_3'),xlabel('Time (sec)')
subplot(4,1,4),hold on,plot(t,qtrue(:,4),'LineWidth',2),plot(t,q3(4,:),'r','LineWidth',1),ylabel('q_4'),xlabel('Time (sec)')

%% Quaternion Error

for i = 1:size(q1,2)
    q1_error(:,i) = qmult(qconj(qtrue(i,:)'))*q1(:,i);
    e1 = Q2phi(q2Q(q1_error(:,i)));
    th_e1(i) = 180/pi*norm(e1);
    th_e1_av(i) = mean(th_e1);
    
    q2_error = qmult(qconj(qtrue(i,:)'))*q2(:,i);
    e2 = unhat(logm(q2Q(q2_error)));
    th_e2(i) = 180/pi*norm(e2);
    th_e2_av(i) = mean(th_e2);
    
    q3_error = qmult(qconj(qtrue(i,:)'))*q3(:,i);
    e3 = unhat(logm(q2Q(q3_error)));
    th_e3(i) = 180/pi*norm(e3);
    th_e3_av(i) = mean(th_e3);
    
end

figure
subplot(3,1,1)
plot(t,th_e1,'.','LineWidth',2)
title('SVD ERRORS')

subplot(3,1,2)
plot(t,th_e2,'.','LineWidth',2)
ylabel('Error (Degrees)')
title('DAV ERRORS')

subplot(3,1,3)
plot(t,th_e3,'.','LineWidth',2)
xlabel('Time (sec)')
title('TRIAD ERRORS')

%% MEKF Application

rN = [rn1, rn2];
yhist = [rb1; rb2];

q0 = q1(:,1);

W = 0.001*eye(6);
V = 0.004*eye(6);

beta = zeros(3,1);
btrue = [];

w_true = 1;

if w_true == 1
    for i = 1:size(y3,1)
        beta = beta + (randn(1,3)*chol(GYRO_ARW^2*eye(3)))';
        whist(:,i) = y3(i,7:9)';
        btrue = [btrue, beta];
    end
else
    for i = 1:size(y3,1)
        beta = beta + (randn(1,3)*chol(GYRO_ARW^2*eye(3)))';
        whist(:,i) = y3(i,7:9)' + beta + (randn(1,3)*chol(GYRO_RND^2*eye(3)))';
        btrue = [btrue, beta];
    end
end

x0 = [q0; 0; 0; 0]; %Initialize with no bias
P0 = (10*pi/180)^2*eye(6); %10 deg. and 10 deg/sec 1-sigma uncertainty

dt = stepSize;

[xhist,Phist] = mekf(x0,P0,W,V,rN,whist,yhist,dt);

qtrue = qtrue';
e = zeros(size(qtrue));
for k = 1:size(qtrue,2)
    e(:,k) = qmult(qconj(qtrue(:,k)))*xhist(1:4,k);
end

%----- Plots -----%
figure;
subplot(4,1,1);
plot(t,qtrue(1,:));
hold on;
plot(t,xhist(1,:));
title('Attitude');
legend('True', 'Estimated');
subplot(4,1,2);
plot(t,qtrue(2,:));
hold on;
plot(t,xhist(2,:));
axis([0 t(end) -1 1])
subplot(4,1,3);
plot(t,qtrue(3,:));
hold on;
plot(t,xhist(3,:));
axis([0 t(end) -1 1])
subplot(4,1,4);
plot(t,qtrue(4,:));
hold on;
plot(t,xhist(4,:));
axis([0 t(end) -1 1])

figure;
subplot(3,1,1);
plot(t,(180/pi)*e(1,:));
hold on
plot(t,(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r','LineWidth',2);
plot(t,-(180/pi)*sqrt(squeeze(Phist(1,1,:))),'r','LineWidth',2);
title('Attitude Error');
axis([0 t(end) -1 1])
subplot(3,1,2);
plot(t,(180/pi)*e(2,:));
hold on
plot(t,(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r','LineWidth',2);
plot(t,-(180/pi)*sqrt(squeeze(Phist(2,2,:))),'r','LineWidth',2);
ylabel('degrees');
axis([0 t(end) -1 1])
subplot(3,1,3);
plot(t,(180/pi)*e(3,:));
hold on
plot(t,(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r','LineWidth',2);
plot(t,-(180/pi)*sqrt(squeeze(Phist(3,3,:))),'r','LineWidth',2);
axis([0 t(end) -1 1])

figure;
subplot(3,1,1);
plot(t,xhist(5,:)-btrue(1,:));
hold on
plot(t,2*sqrt(squeeze(Phist(4,4,:))),'r','LineWidth',2);
plot(t,-2*sqrt(squeeze(Phist(4,4,:))),'r','LineWidth',2);
title('Bias Error');
axis([0 t(end) -1 1])
subplot(3,1,2);
plot(t,xhist(6,:)-btrue(2,:));
hold on
plot(t,2*sqrt(squeeze(Phist(5,5,:))),'r','LineWidth',2);
plot(t,-2*sqrt(squeeze(Phist(5,5,:))),'r','LineWidth',2);
axis([0 t(end) -1 1])
subplot(3,1,3);
plot(t,xhist(7,:)-btrue(3,:));
hold on
plot(t,2*sqrt(squeeze(Phist(6,6,:))),'r','LineWidth',2);
plot(t,-2*sqrt(squeeze(Phist(6,6,:))),'r','LineWidth',2);
axis([0 t(end) -1 1])
% close all;

%% Helper Functions

function drdt = orb_prop(r,mu)
drdt = [r(4:6); -mu*r(1:3)/norm(r(1:3))^3];
end

function drdt = J2_prop(r,mu,J2,Re)

Z = 5*(r(3)^2)/(norm(r(1:3))^2);

drdt = [r(4:6);...
    -mu*r(1)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-1));...
    -mu*r(2)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-1));...
    -mu*r(3)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-3))];

end

function [a, e, inc, RAAN, w, M, E, anom, revs_per_day] =...
    TLE_Reader(filename)
% OE2ECI Converts orbital elements to r, v in inertial frame
%
%   Notes:
%       In cases of equatorial and/or circular orbits, it is assumed
%       that valid orbital elements are provided as inputs (ie. there is
%       no back-end validation)
%
% Inputs:
%      filename - Two Line Element (TLE) .txt file
%
% Outputs:
%      a - semi-major axis of orbit [km]
%      e - eccentricity of orbit
%      inc - inclination of orbit [deg]
%      RAAN - right ascension of the ascending node [deg]
%      w - argument of periapsis [deg]
%      M - mean anomaly [deg]
%      E - eccentric anomaly [deg]
%      anom - true anomaly [deg]
%      revs_per_day - mean motion

fid = fopen(filename, 'rb');
L1c = fscanf(fid,'%21c%',1);
L2c = fscanf(fid,'%71c%',1);
L3c = fscanf(fid,'%71c%',1);
fprintf(L1c);
fprintf(L2c);
fprintf([L3c,'\n']);
fclose(fid);

fid = fopen(filename, 'rb');
L1 = fscanf(fid,'%24c%*s',1);
L2 = fscanf(fid,'%d %6d %*c%5d%*3c%f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
L3 = fscanf(fid,'%d%6d%f%f%f%f%f%f%f',[1,8]);
fclose(fid);

date  = L2(1,4);              % Epoch Date and Julian Date Fraction
Db    = L2(1,5);             % Ballistic Coefficient
inc   = L3(1,3);             % Inclination [deg]
RAAN  = L3(1,4);             % Right Ascension of the Ascending Node [deg]
e     = L3(1,5)/1e7;         % Eccentricity
w     = L3(1,6);             % Argument of periapsis [deg]
M     = L3(1,7);             % Mean anomaly [deg]
n     = L3(1,8);             % Mean motion [Revs per day]

% Orbital elements
mu = 398600; %  Standard gravitational parameter for the earth
revs_per_day = n;
n = revs_per_day*2*pi/(24*3600);
a = (mu/n^2)^(1/3);
E = M2E(M*pi/180,e,10^-6);  %[rad]
anom = E2anom(E, e);        %[rad]
E = E*180/pi;               %[deg]
anom = anom*180/pi;         %[deg]

% Six orbital elements
OE = [a e inc RAAN w M];

%Date Conversion
if (date<57000)
    epoch = datestr(datenum(round(date/1000)+2000,1,0)+mod(date,1000));
else
    epoch = datestr(datenum(round(date/1000)+1900,1,0)+mod(date,1000));
end

fprintf('\n a[km]     e        inc[deg]   RAAN[deg]   w[deg]     M[deg]')
fprintf('\n %4.2f   %4.4f   %4.4f    %4.4f    %4.4f   %4.4f', OE);
fprintf('\n\nEpoch: %s UTC', epoch)
fprintf('\n\n---------- End of TLE Import ----------\n')

end

function [r_eci, v_eci] = OE2ECI(a, e, i, Om, w, anom, mu)
% OE2ECI Converts orbital elements to r, v in inertial frame
%
%   Notes:
%       In cases of equatorial and/or circular orbits, it is assumed
%       that valid orbital elements are provided as inputs (ie. there is
%       no back-end validation)
%
% Inputs:
%      a - semi-major axis of orbit [km]
%      e - eccentricity of orbit
%      i - inclination of orbit [deg]
%      Om - right ascension of the ascending node [deg]
%      w - argument of periapsis [deg]
%      anom - true anomaly [deg]
%      mu - central body gravitational parameters [km^3/s^2]
%
% Outputs:
%   r_eci - 3x1 vector of radius in ECI frame [km]
%   v_eci - 3x1 vector of velocity in ECI frame [km/s]

n = sqrt(mu/a^3);      % rad/s

E = anom2E(deg2rad(anom), e);    % rad

% Compute radius and velocity of orbit in perifocal coordinates
rPeri = [        a*(cos(E) - e);
    a*sqrt(1 - e^2)*sin(E);
    0];

vPeriComp = [             -sin(E);
    sqrt(1 - e^2)*cos(E);
    0];
vPeri = (a*n)/(1 - e*cos(E))*vPeriComp;

% Develop rotation matrix depending on orbit shape/inclination
if i == 0 && e ~= 0         % Equatorial + elliptical
    rotPeri2ECI = rotz(w);
elseif e == 0 && i ~= 0     % Circular + inclined
    rotPeri2ECI = rotz(Om)*rotx(i);
elseif i == 0 && e == 0     % Equatorial + circular
    rotPeri2ECI = 1;
else                        % Elliptical + inclined
    rotPeri2ECI = rotz(Om)*rotx(i)*rotz(w);
end

% Rotate vectors into ECI frame
r_eci = rotPeri2ECI*rPeri;
v_eci = rotPeri2ECI*vPeri;
end

function v = unhat(M)
v(1) = -M(2,3);
v(2) = M(1,3);
v(3) = -M(1,2);
end

function M = hat(v)
M = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
end

function q = Q2quat(Q)
phi = unhat(logm(Q));
theta = norm(phi);
r = (phi/theta)';
q = [r*sin(theta/2); cos(theta/2)];
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

% Matrix Linear Equation
x_dot = [vvec; -mu*rvec/norm(rvec)^3];
x_dot(7:9,1) = -J\cross(om0,J*om0 + rho);
x_dot(10:13,1) = (1/2)*qhat*[om0; 0];

end
