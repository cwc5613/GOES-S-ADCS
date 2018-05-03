clear;
clc;
close all;

[a, e, inc, RAAN, w, M0, E, anom, revs_per_day] =...
    TLE_Reader('SKYSAT1_TLE.txt');

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
tol = 1e-6;

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci; v_eci];

orbitCount = 500;
stopTime = orbitCount/revs_per_day*24*60*60;
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

function dav_q_method(rB,rN,Q_true)
B = zeros(3);
Z = zeros(3,1);
for i=1:size(rB,2)
    B = B + rB(:,i)*rN(:,i)';
    Z = Z + cross(rB(:,i),rN(:,i));
end

K = [B + B' - trace(B)*eye(3) Z;Z' trace(B)];
[eigvec, eigval] = eig(K);
[eigmax, idx] = max(diag(eigval));

q_NB = eigvec(:,idx)/norm(eigvec(:,idx))
v = q_NB(1:3);
s = q_NB(4);

Q_NB_DAVEN = eye(3) + 2*hat(v)*(hat(v) + s*eye(3));
Q_error = Q_true'*Q_NB_DAVEN;
e = unhat(logm(Q_error));
th_e = 180/pi*norm(e);
end

function v = unhat(M)
v(1) = -M(2,3);
v(2) = M(1,3);
v(3) = -M(1,2);
end

function M = hat(v)
M = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
end

