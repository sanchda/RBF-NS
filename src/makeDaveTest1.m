function U = makeDaveTest1(N0, X, t, nu)
% AUTHOR:   David Sanchez
% DATE:     July 2013
% MODIFIED: 7/5/2013

% Returns u(x,t) evaluated at all x in X at the prescribed t, where u is
% the reference solution I made up, based on Ganesh 2009.

% u is defined to be g(t)*core
% g(t) = nu*exp(-t)
% core = sum Y_0^L + 2*Y_{1:L}^L, for L = 1:N0 ; ^L a superscript, not
% exponent

% Notes:
% indexing through the output of getDivFree might seem strange.  Check
% getDivFree.m for details on the output formatting.

%==========================================================================
%                             Inline defines                       
%==========================================================================
g  = @(t) nu*exp(-t);

% Just to zero out U with something!
Z = getDivFree(1,X);
W1 = Z(:,4:6) + 2*Z(:,7:9);

% At first, U holds the running total for the sum Y_... terms.
U = 0*W1;

%==========================================================================
%                            svharms loop                       
%==========================================================================
for L = 1:N0
    Z = getDivFree(L,X);
    Z = Z(:,(3*L+1):end);
    U = U + Z(:,1:3) + ...
        2*[sum(Z(:,4:3:end),2) sum(Z(:,5:3:end),2) sum(Z(:,6:3:end),2)];
end

    U = t*g(t)*U;
end

