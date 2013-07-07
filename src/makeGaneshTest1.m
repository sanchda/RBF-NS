function U = makeGaneshTest1(N0, X, t, nu)
% AUTHOR:   David Sanchez
% DATE:     December 2012
% MODIFIED: 7/5/2013

% Returns u(x,t) evaluated at all x in X at the prescribed t, where u is
% the reference solution in Ganesh 2009.

% u is defined to be t*g(t)*core + g(t)*W1 + (1-t)g(t)*W2
% g(t) = nu*exp(-t)(sin(5t) + cos(10t))
% W1   = Z_0^1 + 2Y_1^1
% W2   = Y_0^2 + 2Y_1^2 + 2Y_2^2
% core = sum Y_0^L + 2*Y_{1:L}^L, for L = 1:N0 ; ^L a superscript, not
% exponent

% Notes:
% indexing through the output of getDivFree might seem strange.  Check
% getDivFree.m for details on the output formatting.

%==========================================================================
%                             Inline defines                       
%==========================================================================
g  = @(t) nu*exp(-t)*(sin(5*t)+cos(10*t));

Z = getDivFree(1,X);
W1 = Z(:,4:6) + 2*Z(:,7:9);

Z = getDivFree(2,X);
W2 = Z(:,7:9) + 2*(Z(:,10:12) + Z(:,13:15));

W2 = 0*W2;

% At first, U holds the running total for the sum Y_... terms.
U = 0*W2;

%==========================================================================
%                            svharms loop                       
%==========================================================================
for L = 1:N0
    Z = getDivFree(L,X);
    Z = Z(:,(3*L+1):end); % need only nonnegative-indexed VHS
    U = U + Z(:,1:3) + ...
        2*[sum(Z(:,4:3:end),2) sum(Z(:,5:3:end),2) sum(Z(:,6:3:end),2)];
end

    U = 0*t*g(t)*U + g(t)*W1 + (t-1)*g(t)*W2;
end

