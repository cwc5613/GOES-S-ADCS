clear;
clc;
close all;

tic

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

stopTime = 2;
stepSize = 1/50;
tspan = [0:stepSize:stopTime];
opts  = odeset('reltol', tol, 'abstol', tol);

[t,y] = ode45(@(t,y) orb_prop(y,mu), tspan, initCond, opts);

%% Quaternions

% Initial Attitude Conditions
Q0 = [0 -1 0;1 0 0;0 0 1];
q0 = Q2quat(Q0);
om0 = [0;1;10.7];
x0 = [r_eci;v_eci;om0;q0];

% ODE45 Solver for Orbital Dynamics
[t, y3] = ode45(@(t,y3)HSTdynamics(mu,J,rho,y3), tspan, x0, opts);
qtrue = y3(:,10:13);

total_sims = 2500

for jj = 1:total_sims
    
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
    
    svd_err_av(jj) = mean(th_e1);
    dav_err_av(jj) = mean(th_e2);
    tri_err_av(jj) = mean(th_e3);
    
end

toc

tri_mean = mean(tri_err_av)
svd_mean = mean(svd_err_av)
dav_mean = mean(dav_err_av)

figure
subplot(3,1,1),plot(tri_err_av,'.')
subplot(3,1,2),plot(svd_err_av,'.')
subplot(3,1,3),plot(dav_err_av,'.')

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