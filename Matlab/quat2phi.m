function phi = quat2phi(q)
q(isnan(q)) = 0;
Q = quat2Q(q);
phi = unhat(logm(Q));
end