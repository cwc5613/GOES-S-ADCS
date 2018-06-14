function [mu, sigma] = mekf(x0, P0, Q, R, rN, whist, yhist, dt)

mu = zeros(7,size(yhist,2));
mu(:,1) = x0;

sigma = zeros(6,6,size(yhist,2));
sigma(:,:,1) = P0;

for k = 1:(size(yhist,2)-1) 
    [mu_pred, A] = prediction(mu(:,k),whist(:,k),dt);
    sig_pred = A*sigma(:,:,k)*A' + Q;   
    [yp, C] = measurement(mu_pred(1:4),rN);   
    
    % Innovation
    inno_quat = quat2phi(qmult(qconj(mu_pred(1:4)))*yhist(1:4,k+1)); 
    inno_body = yhist(5:end,k+1) - yp(5:end); 
    inno = [inno_quat; inno_body];
    S = C*sig_pred*C' + R; 

    % Kalman Gain
    K = sig_pred*C'/S; 
    
    % Update
    dx = K*inno;
    phi = dx(1:3);
    dq = phi2quat(phi);
    dq = dq/norm(dq);
    
    % Quaternion and Bias Update
    mu(1:4,k+1) = qmult(mu_pred(1:4))*dq;
    mu(5:7,k+1) = mu_pred(5:7) + dx(4:6);
    
    % Covariance Update
    sigma(:,:,k+1) = (eye(6) - K*C)*sig_pred*(eye(6) - K*C)' + K*R*K';
end
end

