function [xhist, Phist] = mekf(x0, P0, W, V, rN, whist, yhist, dt)

xhist = zeros(7,size(yhist,2));
xhist(:,1) = x0;

Phist = zeros(6,6,size(yhist,2));
Phist(:,:,1) = P0;


for k = 1:(size(yhist,2)-1) 
    [x_p, A] = prediction(xhist(:,k),whist(:,k),dt);
    P_p = A*Phist(:,:,k)*A' + W;   
    [yp, C] = measurement(x_p(1:4),rN);   
    
    % Innovation
    z = yhist(:,k) - yp;
    S = C*P_p*C' + V;

    % Kalman Gain
    L = P_p*C'*S^-1;
    
    % Update
    dx = L*z; % [dphi; dbeta], 6x1
    phi = dx(1:3);
    dq = [1/2*phi ; 1 - 1/8*phi'*phi];
    dq = dq/norm(dq);
    
    % Quaternion and Bias Update
    xhist(1:4,k+1) = qmult(x_p(1:4))*dq;
    xhist(5:7,k+1) = x_p(5:7) + dx(4:6);
    
    % Covariance Update
    Phist(:,:,k+1) = (eye(6) - L*C)*P_p*(eye(6) - L*C)' + L*V*L';

end
end

