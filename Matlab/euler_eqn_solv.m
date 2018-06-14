function wdot = euler_eqn_solv(t,w,J)
wdot = J\-(cross(w,(J*w)));
end