clear;
clc;
close all;

TLE_Reader;

mu = 398600;
J2 = 1.081874*10^-3;
Re = 6378137*10^-3;
M0 = M;
tol = 10^-5;

[r_eci, v_eci] = OE2ECI(a, e, inc, RAAN, w, anom, mu);

initCond = [r_eci, v_eci];

f0 = 0;
v0 = v_eci;
r0 = r_eci;

options = simset('SrcWorkspace','current');
orbitCount = 500;
stopTime = orbitCount/revs_per_day*24*60*60;
stepSize = 100;
sim('NUM_PROP',[],options);

out_r_NUM_PROP = reshape(NUM_r,[size(NUM_r,1),size(NUM_r,3)]);
out_v_NUM_PROP = reshape(NUM_v,[size(NUM_v,1),size(NUM_v,3)]);

figure
hold on
plot3(out_r_NUM_PROP(1,:),out_r_NUM_PROP(2,:),out_r_NUM_PROP(3,:),...
    'c','LineWidth',1)
earthPlot(1)
axis equal
hold off
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

sim('J2_Accel',[],options);

out_r_J2_Accel = reshape(J2_r,[size(J2_r,1),size(J2_r,3)]);
out_v_J2_Accel = reshape(J2_v,[size(J2_v,1),size(J2_v,3)]);

figure
hold on
plot3(out_r_J2_Accel(1,:),out_r_J2_Accel(2,:),out_r_J2_Accel(3,:),...
    'c','LineWidth',1)
earthPlot(1)
axis equal
hold off
xlabel('[km]')
ylabel('[km]')
zlabel('[km]')

