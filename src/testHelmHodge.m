function [ Maxerr, L2err, kappa ] = testHelmHodge(x, W, U0div, U0crl, H, epsilon)
% Tests the divergence-free code.

    N = size(U0div,1);

    % Combine the input vector fields
    U0 = U0div + U0crl;
    
    % Zonal and meridional bases
    d = @(x) (1/sqrt(1-x(3)^2))*[(-x(3)*x(1)) (-x(3)*x(2)) (1-x(3)^2)];
    e = @(x) (1/sqrt(1-x(3)^2))*[(-x(2)) x(1) 0];
    
    % The function Q(x) generates the matrix which will project a vector in R3
    % onto the tangent space of the sphere at x \in S2
    Q = @(x) [0 x(3) (-x(2)); (-x(3)) 0 x(1);  x(2) (-x(1)) 0];
    
    % P(x) projects into the normal space
    P = @(x) eye(3) - [x(1);x(2);x(3)]*[x(1) x(2) x(3)];
   
    % Construct the SBF thingy
    PSIdiv  = @(x,y) (Q(x)')*(-H(x,y,epsilon))*Q(y);
    PSIcrl  = @(x,y) (P(x)')*(-H(x,y,epsilon))*P(y);
    PSI     = @(x,y) PSIdiv(x,y) + PSIcrl(x,y);
    
    Adiv =  makeAdiv(x, PSI, d, e);
    
    kappa = cond(Adiv);

    % TODO: build gamdel in a vectorized way.
    gamdel = zeros(2*N,1);
    buff = zeros(2,3);

    for i=1:N
        buff(1,:) = d(x(i,:))';
        buff(2,:) = e(x(i,:))';
        gamdel(2*i-1:2*i) = buff*U0(i,:)';
    end
    
    albet = Adiv\gamdel;
    
    dmat = zeros(N,3);
    emat = zeros(N,3);
    
    for i = 1:N
        dmat(i,:) = d(x(:,1));
        emat(i,:) = e(x(:,1));
    end
   
 
    % recover the interpolation coefficients
    interpcoeffs = zeros(N,3);
    for i = 1:N
        interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
    end
    
    % TODO:  improve performance
    Udivfree = zeros(N,3);
    Ucrlfree = zeros(N,3);
    for i = 1:N
        for j = 1:N
            Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
            Ucrlfree(i,:) = Ucrlfree(i,:) + (PSIcrl(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
        end
    end
    
    % Figure out the L-infty residual
    err = abs(U0div - Udivfree);
    errsqrt = sqrt(err(:,1).^2 + err(:,2).^2 + err(:,3).^2);
    Maxerr = max(errsqrt);
    
    % TODO:
    % fit epsilon
    % vector lap use div-free vector curl of SVH
    
    % Figure out the L-2 residual
    L2err = 0;
    for i = 1:N
        L2err = L2err + W(i)*(err(i,:)*err(i,:)');
    end
    L2err = sqrt(L2err);
    
end

