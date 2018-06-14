function [K,S] = tvlqr(A,B,Q,R,S)
K = (R+B'*S*B)\B'*S*A;
S = Q+K'*R*K+(A-B*K)'*S*(A-B*K);