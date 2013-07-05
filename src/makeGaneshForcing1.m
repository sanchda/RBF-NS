function f = makeGaneshForcing1(N0, X, t, nu, grad, Pxmat)
% Returns f(x,t) evaluated at all x in X at the prescribed t, where f is
% the forcing corresponding to the reference solution in Ganesh 2009.

% u is defined to be t*g(t)*core + g(t)*W1 + (1-t)g(t)*W2
% g(t) = nu*exp(-t)(sin(5t) + cos(10t))
% W1   = Z_0^1 + 2Y_1^1
% W2   = Y_0^2 + 2Y_1^2 + 2Y_2^2
% core = sum Y_0^L + 2*Y_{1:L}^L, for L = 1:N0

% f = u_t - nu*lapu + covdu

%==========================================================================
%                             Inline defines                       
%==========================================================================
g  = @(t)  nu*exp(-t)*(sin(5*t)+cos(10*t));
gp = @(t)  nu*exp(-t)*(-sin(5*t)-cos(10*t) +5*cos(5*t) -10*sin(10*t));

Z = getDivFree(1,X);
W1 = [Z(:,4) Z(:,5) Z(:,6)] + 2*[Z(:,7) Z(:,8) Z(:,9)];

Z = getDivFree(2,X);
W2 = [Z(:,7) Z(:,8) Z(:,9)] + 2*[Z(:,10) Z(:,11) Z(:,12)] + 2*[Z(:,13) Z(:,14) Z(:,15)];

U = 0*X;
Ulap = 0*U;

%==========================================================================
%                            svharms loop                       
%==========================================================================
for L = 1:N0
    Z = getDivFree(L,X);
    Z = Z(:,(3*L+1):end);
    Ubuff = Z(:,1:3) + [2*sum(Z(:,4:3:end),2),2*sum(Z(:,5:3:end),2),2*sum(Z(:,6:3:end),2)];
    Ulap = Ulap - L*(L+1)*Ubuff;
    U = U + Ubuff;
end

%==========================================================================
%                             Make Forcing                       
%==========================================================================
    Ulap = t*g(t)*Ulap - g(t)*1*(1+1)*W1 - (t-1)*g(t)*2*(2+1)*W2;
    Ut   = (g(t) + t*gp(t))*U + gp(t)*W1 + ((t-1)*gp(t) - g(t))*W2;
    U = t*g(t)*U + g(t)*W1 + (t-1)*g(t)*W2;
    
    % Compute the covariant derivative
    % TODO: write this out directly in terms of the spherical harmonics
    covu = grad*U(:,1);
    covu = reshape(covu,[],3);
    covu = U(:,1).*covu(:,1) + U(:,2).*covu(:,2) + U(:,3).*covu(:,3);

    covv = grad*U(:,2);
    covv = reshape(covv,[],3);
    covv = U(:,1).*covv(:,1) + U(:,2).*covv(:,2) + U(:,3).*covv(:,3);

    covw = grad*U(:,3);
    covw = reshape(covw,[],3);
    covw = U(:,1).*covw(:,1) + U(:,2).*covw(:,2) + U(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    f = Ut + covU - nu*Ulap;
end

