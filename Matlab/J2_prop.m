function drdt = J2_prop(r,mu,J2,Re)

Z = 5*(r(3)^2)/(norm(r(1:3))^2);

drdt = [r(4:6);...
    -mu*r(1)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-1));...
    -mu*r(2)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-1));...
    -mu*r(3)/(norm(r(1:3))^3)*(1-3/2*J2*(Re/norm(r(1:3)))^2*(Z-3))];

end