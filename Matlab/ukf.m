function [mu, sigma] = ukf(x0, P0, Q, R, rN, whist, yhist, dt)

mu = zeros(7,size(yhist,2));
mu(:,1) = x0;

sigma = zeros(6,6,size(yhist,2));
sigma(:,:,1) = P0;

for k = 1:(size(yhist,2)-1)
    
    %--------------------------------------------
    % Calculate Sigma Points
    %--------------------------------------------
    
    [X,wm,wc] = sigpoint_calc(mu(:,k),sigma(:,:,k));
    
    %--------------------------------------------
    % Prediction Phase
    %--------------------------------------------
    
    for i = 1:size(X,2)
        [Xstar(:,i),~] = prediction(X(:,i),whist(:,k),dt);
    end
    
    [mu_pred,W_prime] = wmean(Xstar,wm);
    Pk = W_prime*diag(wc)*W_prime' + Q;
    
    %--------------------------------------------
    % Update Phase
    %--------------------------------------------
    
    [X,wm,wc] = sigpoint_calc(mu_pred,Pk);
    
    for i = 1:size(X,2)
        [yp(:,i), ~] = measurement(X(1:4,i),rN);
    end
    
    [y_est,W_primeY] = wmean(yp,wm);
    [~,W_primeX] = wmean(X,wm);
        
    Pzz = W_primeY*diag(wc)*W_primeY' + R;
    Pxz = W_primeX*diag(wc)*W_primeY';

    K = Pxz/Pzz;
    
    % Innovation
    inno_quat = quat2phi(qmult(qconj(y_est(1:4)))*yhist(1:4,k+1)); 
    inno_body = yhist(5:end,k+1) - y_est(5:end); 
    inno = [inno_quat; inno_body];
    
    % Update
    dx = K*inno;
    phi = dx(1:3);
    dq = phi2quat(phi);
    dq = dq/norm(dq);
    
    % Quaternion and Bias Update
    mu(1:4,k+1) = qmult(mu_pred(1:4))*dq;
    mu(5:7,k+1) = mu_pred(5:7) + dx(4:6);
    
    % Covariance Update
    sigma(:,:,k+1) = (Pk - Pxz/Pzz*Pxz');
end
end

function [X,wm,wc] = sigpoint_calc(mu,sig)

L = length(mu);
lam = 2;

gamma = sqrt(L+lam);

sig(isnan(sig)) = 0;
S = gamma*chol(sig,'L');

n = length(S);
q = mu(1:4);
w = mu(5:end);

Xq(:,1) = mu(1:4);
Xw(:,1) = mu(5:7);
wm(1) = lam/(L+lam);
wc(1) = lam/(L+lam) + (1-(1e-3)^2+2);

for j = 2:n+1
    W = S;
    ew = W(1:3,j-1);
    qw = phi2quat(ew);
    Xq(:,j) = qmult(q)*qw;
    wm(j) = 1/(2*(L+lam));
    wc(j) = 1/(2*(L+lam));
end

for j = n+2:(2*n)+1
    W = -S;
    ew = W(1:3,j-(n+1));
    qw = phi2quat(ew);
    Xq(:,j) = qmult(q)*qw;
    wm(j) = 1/(2*(L+lam));
    wc(j) = 1/(2*(L+lam));
end

for j = 2:n+1
    W = S;
    Xw(:,j) = w+W(4:6,j-1);
end

for j = n+2:(2*n)+1
    W = -S;
    Xw(:,j) = w+W(4:6,j-(n+1));
end

X = [Xq;Xw];
wm = wm/sum(wm);
wc = wc/sum(wc);

end
function [pred,inno] = wmean(sigpts,wm)

M = zeros(4);

for i = 1:size(sigpts,2)
    xhat2(:,i) = wm(:,i)*sigpts(5:end,i);
    M = M + wm(i)*sigpts(1:4,i)*sigpts(1:4,i)';
end
K = 4*M-sum(wm)*eye(4);
[v,d] = eig(K);
[~,I] = max(diag(d));
mean_q = v(:,I)/norm(v(:,I));

if dot(mean_q,sigpts(1:4,1)) < 0
    mean_q = -mean_q;
end

mean_2 = sum(xhat2,2);
pred = [mean_q;mean_2];

for i = 1:size(sigpts,2)
    wW(:,i) = (sigpts(5:end,i) - pred(5:end));
    rW(:,i) = quat2phi(qmult(qconj(pred(1:4)))*sigpts(1:4,i));
end

inno = [rW;wW];

end