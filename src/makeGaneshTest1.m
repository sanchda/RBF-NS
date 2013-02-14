function U = makeGaneshTest1(N0, X, t, nu)
% Returns u(x,t) evaluated at all x in X at the prescribed t, where u is
% the reference solution in Ganesh 2009.

% u is defined to be t*g(t)*core + g(t)*W1 + (1-t)g(t)*W2
% g(t) = nu*exp(-t)(sin(5t) + cos(10t))
% W1   = Z_0^1 + 2Y_1^1
% W2   = Y_0^2 + 2Y_1^2 + 2Y_2^2
% core = sum Y_0^L + 2*Y_{1:L}^L, for L = 1:N0

%==========================================================================
%                          Inline defines                       
%==========================================================================
g  = @(t) nu*exp(-t)*(sin(5*t)+cos(10*t));

[Z1x, Z1y, Z1z] = getDivFree(1,X);

W1 = dsph(1,X(:,1),X(:,2),X(:,3)); % gather up the spharms
W1 = W1(:,2) + 2*W1(:,3);          % add them together
W2 = dsph(2,X(:,1),X(:,2),X(:,3));
W2 = W2(:,3) + 2*sum(W2(:,4:5),2);

U = 0*X;


%==========================================================================
%                          Spharm loop                       
%==========================================================================
for L = 1:N0
    s2Harms = dsph(L,X(:,1),X(:,2),X(:,3));
    U = U + sum(s2Harms(:,(L+2):end),2) + s2Harms(:,L+1);
end

    U = t*g(t)*core + g(t)*W1 + (t-1)*g(t)*W2;
end

