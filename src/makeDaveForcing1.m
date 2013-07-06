function f = makeDaveForcing1(N0, X, t, nu, projgrad, Pxmat)
% AUTHOR:   David Sanchez
% DATE:     July 2013
% MODIFIED: 7/5/2013

% Returns f(x,t) evaluated at all x in X at the prescribed t, where f is
% the forcing corresponding to the reference solution I made up myself, but
% it's based on Ganesh 2009.

% u is defined to be g(t)*core
% g(t) = nu*exp(-t)
% core = sum(Y_0^L + sum(2*Y_{1:L}^L,2), for L = 1:N0)

% f = u_t - nu*lapu + covdu

% Notes:
% indexing through the output of getDivFree might seem strange.  Check
% getDivFree.m for details on the output.

%==========================================================================
%                              Inline defines
%==========================================================================
% g is given by Ganesh; gp is its time-derivative.
% TODO explore simplification of gp to incur fewer operations or evaluate
% fewer trigonometric functions
g  = @(t)  nu*exp(-t);
gp = @(t)  -nu*exp(-t);

% Just to zero shit out
Z = getDivFree(1,X);
W1 = Z(:,4:6) + 2*Z(:,7:9);

U = 0*W1;
Ulap = 0*U;

%==========================================================================
%                      div-free VSH sum-building loop
%==========================================================================
for L = 1:N0
    Z = getDivFree(L,X);
    Z = Z(:,(3*L+1):end); % need only nonnegative-indexed VHS
    Ubuff = Z(:,1:3) + ...
        2*[sum(Z(:,4:3:end),2) sum(Z(:,5:3:end),2) sum(Z(:,6:3:end),2)];
    Ulap = Ulap - (L*(L+1))*Ubuff;  % running total of the summation term for
                                  % synthesis with vector Laplacian
    U = U + Ubuff;                % running total of the summation term for
                                  % normal evaluation
end

%==========================================================================
%                             Make Forcing                       
%==========================================================================
% These quantities were determined by applying the differential operators
% to the analytic solution in Ganesh.  In particular, note that the
% eigenvalues under the Laplacian for Z_{L,M} are -L(L+1), as the Z are
% just div-free vector spherical harmonics.
%
%

    % Application of the vector Laplacian to Ganesh's reference sol'n
    Ulap = t*g(t)*Ulap;
    
    % Application of d/dt
    Ut   = (g(t) + t*gp(t))*U;
    
    % reference solution, for numerical evaluation of covariant derivative
    U = t*g(t)*U;
    
    % Apply covariant Derivative
    % cov_u(u) = Px*[U .* grad(U(:,1); U .* grad(U(:,2); U .* grad(U(:,3)]
    % TODO: write this out directly in terms of the spherical harmonics
    covu = projgrad*U(:,1);
    covu = reshape(covu,[],3);
    covu = U(:,1).*covu(:,1) + U(:,2).*covu(:,2) + U(:,3).*covu(:,3);

    covv = projgrad*U(:,2);
    covv = reshape(covv,[],3);
    covv = U(:,1).*covv(:,1) + U(:,2).*covv(:,2) + U(:,3).*covv(:,3);

    covw = projgrad*U(:,3);
    covw = reshape(covw,[],3);
    covw = U(:,1).*covw(:,1) + U(:,2).*covw(:,2) + U(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    
    % Define the forcing
    f = Ut + covU - nu*Ulap;

end

