function drdt = orb_prop(r,mu)

drdt = [r(4:6); -mu*r(1:3)/norm(r(1:3))^3];

end