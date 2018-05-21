function [q,Q_NB_DAVEN] = davenport_method(rN1,rN2,rN3,rB1,rB2,rB3,v1,v2,v3)

B = zeros(3);
Z = zeros(3,1);

B = B + rB1*rN1'/v1;
B = B + rB2*rN2'/v2;
B = B + rB3*rN3'/v3;

Z = Z + cross(rB1,rN1)/v1;
Z = Z + cross(rB2,rN2)/v2;
Z = Z + cross(rB3,rN3)/v3;

K = [B + B' - trace(B)*eye(3) Z;Z' trace(B)];
[eigvec, eigval] = eig(K);
[~, idx] = max(diag(eigval));

q_NB = eigvec(:,idx)/norm(eigvec(:,idx));
v = q_NB(1:3);
s = q_NB(4);

Q_NB_DAVEN = eye(3) + 2*hat(v)*(hat(v) + s*eye(3)); 
q = Qtoq(Q_NB_DAVEN);

end