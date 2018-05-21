function [q,Q_NB_SVD] = svd_method(rN1,rN2,rN3,rB1,rB2,rB3,v1,v2,v3)

B = zeros(3);

B = B + rB1*rN1'/v1;
B = B + rB2*rN2'/v2;
B = B + rB3*rN3'/v3;

[U,~,V] = svd(B');
Q_NB_SVD = U*[1,0,0;0,1,0;0,0,det(U)*det(V)]*V';
q = Qtoq(Q_NB_SVD);

end

