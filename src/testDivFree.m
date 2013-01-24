function [ Maxerr, L2err, kappa] = testDivFree(x, W, U0, H, eps )
% Tests the divergence-free code.

    N = size(U0,1);

    % Zonal and meridional bases
    d = @(x) (1/sqrt(1-x(3)^2))*[(-x(3)*x(1)) (-x(3)*x(2)) (1-x(3)^2)];
    e = @(x) (1/sqrt(1-x(3)^2))*[(-x(2)) x(1) 0];
    
    % The function Q(x) generates the matrix which will project a vector in R3
    % onto the tangent space of the sphere at x \in S2
    Q = @(x) [0 x(3) (-x(2)); (-x(3)) 0 x(1);  x(2) (-x(1)) 0];
   
    PSI  = @(x,y) Q(x)*(-H(x,y,eps))*(Q(y)');
    
    Adiv =  makeAdiv(x, PSI, d, e);
    
    kappa = cond(Adiv);

    % TODO: build gamdel in a better way.  Also, seriously--did you wake
    % up on the retarded side of the bed when you named this array?
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
   
 
    interpcoeffs = zeros(N,3);
    
    for i = 1:N
        interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
    end
    
    % TODO:  improve performance
    Udivfree = zeros(N,3);
    for i = 1:N
        for j = 1:N
            Udivfree(i,:) = Udivfree(i,:) + (PSI(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
        end
    end
    
    % Figure out the L-infty residual
    err = abs(U0 - Udivfree);
    errsqrt = sqrt(err(:,1).^2 + err(:,2).^2 + err(:,3).^2);
    Maxerr = max(errsqrt);
    
    % Figure out the L-2 residual
    L2err = 0;
    for i = 1:N
        L2err = L2err + W(i)*(err(i,:)*err(i,:)');
    end
    L2err = sqrt(L2err);
    
end

