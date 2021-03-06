function [U,t] = navierstokes(x, U0, h, t, M, nu, omega, N0, lap, projgrad, Aleray, Pxmat)
% AUTHOR:   David Sanchez
% DATE:     August 2012
% MODIFIED: 7/11/2013
    
%==========================================================================
%                               Description
%==========================================================================
% TODO


    
%==========================================================================
%                         Assumptions & Operation
%==========================================================================
%
%==========================WARNING WARNING WARNING=========================
% THIS IS ALL A MASSIVE TODO.  NONE OF THE DEBUG ROUTINES ARE CURRENTLY
% IMPLEMENTED (2/18/2013)
%
% The values in x are assumed to be on the unit sphere, with coordinates
% given in Cartesian x,y,z (i.e., x^2 + y^2 + z^2 = 1).  This code does not
% check where these values lie, but will probably give incorrect solutions
% if this assumption is broken.
%
% U0 is assumed to be a velocity field with values corresponding to the
% values in x.
%
% H is a handle to a function evaluating the Hessian of the desired RBF at
% a point (i.e., H([x1,x2,x3],[y1,y2,y3],epsilon) = [z1,z2,z3]).  This is a
% rather awkward interface and should be able to be fixed by utilizing
% Matlab's symbolic toolkit to compute derivatives.
%
% epsilon is the shape parameter.  Care must be taken to choose the shape
% parameter appropriately for a given size(x,1).
%
% nu, omega are viscocity and coriolis force coefficients, respectively.
% It should be noted that values of nu can dramatically rescale the
% dominance of certain terms (see referenced paper (TODO: write paper))
%
% debugvec is an optional parameter which engages unit tests.  These tests
% are designed to check various properties of the inputs and this
% function's internal routines.
%
% debugvec = [isDivFree, isDFVSH, getLaplacianErr, getHodgeErr]
%
% isDivFree - engages code to verify that the input U0 is divergence-free.
%             will force function to return with as little code execution
%             as possible.  Blocks Navier-Stokes simulation.
%
% isDFVSH   - verifies that the input is a vector spherical harmonic.  will
%             force function to return while executing as little unrelated
%             code as possible.  Blocks Navier-Stokes simulation.
%
% getLaplacianErr - Assumes that the input U0 is a divergence-free vector
%                   spherical harmonic, returning the error in the surface
%                   Laplacian.  Blocks Navier-Stokes simulation.
%
% getHodgeErr     - Checks the error of splitting U0 into curl-free and
%                   divergence-free components.  TODO:  check the recovery
%                   of strictly curl- or div-free vector fields.
%
%


%==========================================================================
%                               Initialization
%==========================================================================
    % Sense the size of the problem from the initial conditions.  The
    % scalars in U0 are assumed to identify with the R3 vectors in X.  If
    % these are not the same length, function should return error.
    % TODO:  return error
%    N = size(U0,1);
    
    % Needed to compute coriolis force
    X = diag(x(:,1));
    Y = diag(x(:,2));
    Z = diag(x(:,3));
    
    % Shortcut for Coriolis force
    zsqrt = 1 - x(:,3).*x(:,3);
    zsqrt = sqrt(zsqrt);

    
%==========================================================================
%                             Begin Timestepping
%==========================================================================  
 
% pre-defines
U = U0;
    
for c = 1:M

    %==============================Begin RK4===============================
    %-----------------------------RK4 Stage 1------------------------------
    arg = U;
        
    % Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    % Covariant Derivative
    % cov_u(u) = Px*[U .* grad(U(:,1); U .* grad(U(:,2); U .* grad(U(:,3)]
    covu = projgrad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = projgrad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = projgrad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    % Coriolis force
    coriolis = [(-Z*arg(:,2) + Y*arg(:,3)) (Z*arg(:,1) - X*arg(:,3)) (-Y*arg(:,1) + X*arg(:,2))];
    coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
    
    % Ganesh force
    fganesh = makeGaneshForcing1(N0, x, t, nu, projgrad, Pxmat);
    
    % Stick it all together
    RK1 = nu*lapU - covU - coriolis + fganesh;
    RK1 = projectDivFree(RK1, Aleray);

    %-----------------------------RK4 Stage 2------------------------------
    arg = U + (h/2)*RK1;
    
    % Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    % Covariant Derivative
    % cov_u(u) = Px*[U .* grad(U(:,1); U .* grad(U(:,2); U .* grad(U(:,3)]
    covu = projgrad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = projgrad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = projgrad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    % Coriolis force
    coriolis = [(-Z*arg(:,2) + Y*arg(:,3)) (Z*arg(:,1) - X*arg(:,3)) (-Y*arg(:,1) + X*arg(:,2))];
    coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;

    % Ganesh force
    fganesh = makeGaneshForcing1(N0, x, t+h/2, nu, projgrad, Pxmat);
    
    %Stick it all together
    RK2 = nu*lapU - covU - coriolis + fganesh;
    RK2 = projectDivFree(RK2, Aleray);
    
    %-----------------------------RK4 Stage 3------------------------------
    arg = U + (h/2)*RK2;
    
    % Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    % Covariant Derivative
    % cov_u(u) = Px*[U .* grad(U(:,1); U .* grad(U(:,2); U .* grad(U(:,3)]
    covu = projgrad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = projgrad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = projgrad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';

    % Coriolis force
    coriolis = [(-Z*arg(:,2) + Y*arg(:,3)) (Z*arg(:,1) - X*arg(:,3)) (-Y*arg(:,1) + X*arg(:,2))];
    coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;  
   
    % Ganesh force
    fganesh = makeGaneshForcing1(N0, x, t+h/2, nu, projgrad, Pxmat);
    
    %Stick it all together
    RK3 = nu*lapU - covU - coriolis + fganesh;
    RK3 = projectDivFree(RK3, Aleray);
    
    %-----------------------------RK4 Stage 4------------------------------
    arg = U + h*RK3;
    
    % Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    % Covariant Derivative
    % cov_u(u) = Px*[U .* grad(U(:,1); U .* grad(U(:,2); U .* grad(U(:,3)]
    covu = projgrad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = projgrad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = projgrad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    % Coriolis force
    coriolis = [(-Z*arg(:,2) + Y*arg(:,3)) (Z*arg(:,1) - X*arg(:,3)) (-Y*arg(:,1) + X*arg(:,2))];
    coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
    
    % Ganesh force
    fganesh = makeGaneshForcing1(N0, x, t+h, nu, projgrad, Pxmat);
    
    %Stick it all together
    RK4 = nu*lapU - covU - coriolis + fganesh;
    RK4 = projectDivFree(RK4, Aleray);

    %============================Stitch Together===========================
    U = U + (h/6)*(RK1 + 2*RK2 + 2*RK3 + RK4);
    
    
    %==========================Determine error=============================
    
    
%     Uganesh = makeGaneshTest1(N0, X, t, nu);
%     
%     error = (Uganesh-U).^2;
%     error = sum(sqrt(sum(error,2)));
%     disp(error);
%     U = sqrt(sum((U - Uganesh).^2,2));

    %==========================Finishing Stuff=============================
    t = t + h;
end

end
