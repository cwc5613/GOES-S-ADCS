function wdot = euler_rotor_solv(t,w,J,rho)
wdot = J\-(cross(w,(J*w+rho)));
end