function U = navierstokes(x, U0, h, M, epsilon, nu, omega, lap, grad, surfeps, Lx, Ly, Lz, Afull, Acrl, PSIfull, PSIcrl, PSIdiv, Pxmat)
%==========================================================================
%                               Description
%==========================================================================
% TODO


    
%==========================================================================
%                         Assumptions & Operation
%==========================================================================
%
%==========================WARNING WARNING WARNING=========================
% THIS IS ALL A MASSIVE TODO.  NONE OF THIS IS CURRENTLY IMPLEMENTED.
% (2/18/2013)
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
    N = size(U0,1);
    
    
    timesMatVec = @(A,b) bsxfun(@times,A,b(:));

    % Zonal and meridional bases.  Works on column-vectors as well as
    % arrays.
    d = @(x) timesMatVec([(-x(:,3).*x(:,1)) (-x(:,3).*x(:,2)) ...
        (1-x(:,3).^2)],(1./sqrt(1-x(:,3).^2)));
    
    e = @(x) timesMatVec([(-x(:,2)) x(:,1) ...
        0*x(:,1)],(1./sqrt(1-x(:,3).^2)));

    
%==========================================================================
%                             Begin Timestepping
%==========================================================================  
 
% pre-defines
U = U0;
gamdel = zeros(2*N,1);
dmat = d(x);
emat = e(x);
t=0;
    
for c = 1:M

    %==============================Begin RK4===============================
    t = t + h;

    %================================RK4 Stage 1===========================
    arg = U;
    
    %Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    %Covariant Derivative
    covu = grad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = grad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = grad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';

    %Stick it all together
    RK1 = nu*lapU - covU;
    RK1 = getDivFree(RK1, dmat, emat, Afull, PSIdiv);

    %================================RK4 Stage 2===========================
    arg = U + 0.5*h*RK1;
    
    %Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    %Covariant Derivative
    covu = grad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = grad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = grad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';

    %Stick it all together
    RK2 = nu*lapU - covU;
    RK2 = getDivFree(RK2, dmat, emat, Afull, PSIdiv);
    
    %================================RK4 Stage 3===========================
    arg = U + 0.5*h*RK2;
    %Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    %Covariant Derivative
    covu = grad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = grad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = grad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';

    %Stick it all together
    RK3 = nu*lapU - covU;
    RK3 = getDivFree(RK3, dmat, emat, Afull, PSIdiv);
    
    %================================RK4 Stage 4===========================
    arg = U + h*RK3;
    
    %Vector Laplacian
    lapU = lap*reshape(arg,[],1);
    lapU = reshape(lapU,[],3);

    %Covariant Derivative
    covu = grad*arg(:,1);
    covu = reshape(covu,[],3);
    covu = arg(:,1).*covu(:,1) + arg(:,2).*covu(:,2) + arg(:,3).*covu(:,3);

    covv = grad*arg(:,2);
    covv = reshape(covv,[],3);
    covv = arg(:,1).*covv(:,1) + arg(:,2).*covv(:,2) + arg(:,3).*covv(:,3);

    covw = grad*arg(:,3);
    covw = reshape(covw,[],3);
    covw = arg(:,1).*covw(:,1) + arg(:,2).*covw(:,2) + arg(:,3).*covw(:,3);

    % Pxmat acts on the row-vectorized form, so the transposition below is
    % necessary.
    covU = Pxmat*reshape([covu covv covw]',[],1);
    covU = reshape(covU,3,[])';
    
    %Stick it all together
    RK4 = nu*lapU - covU;
    RK4 = getDivFree(RK4, dmat, emat, Afull, PSIdiv);

    %============================Stitch Together===========================

    U = U + (1/6)*(RK1 + 2*RK2 + 2*RK3 + RK4);

%U = Uthis;
    
%     
%     %=============================RK4 part 1===============================
%     % Compute the Laplacian
%     lapu = -[(lapux*U(:,1)+lapuy*U(:,2)+lapuz*U(:,3)) (lapvx*U(:,1)+lapvy*U(:,2)+lapvz*U(:,3)) (lapwx*U(:,1)+lapwy*U(:,2)+lapwz*U(:,3))];
% 
%     % Compute the covariant derivative
%     covU = -[(covdux*U(:,1)+covduy*U(:,2)+covduz*U(:,3)) (covdvx*U(:,1)+covdvy*U(:,2)+covdvz*U(:,3)) (covdwx*U(:,1)+covdwy*U(:,2)+covdwz*U(:,3))];
%     covrep = (U(:,1).*(Lx*U(:,1)) + U(:,2).*(Ly*U(:,2)) + U(:,3).*(Lz*U(:,3)));
%     covU = covU + repmat(covrep,1,3);
%     
%     % Compute the coriolis force
%     coriolis = [(-Z*U(:,2) + Y*U(:,3)) (Z*U(:,1) - X*U(:,3)) (-Y*U(:,1) + X*U(:,2))];
%     coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
%     
%     fU = makeGaneshForcing1(N0, x, t, nu, Lx, Ly, Lz, covdux, covduy, covduz, covdvx, covdvy, covdvz, covdwx, covdwy, covdwz);
%     
%     U = -covU + nu*lapu - coriolis + fU;
% U = -covU + nu*lapu;
%     %-----Leray Projection
%     % TODO: build gamdel in a vectorized way.
% 
%     for i=1:N
%         buff(1,:) = d(x(i,:))';
%         buff(2,:) = e(x(i,:))';
%         gamdel(2*i-1:2*i) = buff*U(i,:)';
%     end
%     
%     albet = Adiv\gamdel;
%    
%     % recover the interpolation coefficients
%     interpcoeffs = zeros(N,3);
%     for i = 1:N
%         interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
%     end
%     
%     % TODO:  improve performance
%     Udivfree = zeros(N,3);
%     for i = 1:N
%         for j = 1:N
%             Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
%         end
%     end
%     
%     U = Udivfree;
%     k1 = U;
%     
%     U = Ubuff + (h/2)*k1;
%     %=============================RK4 part 2===============================
%     % Compute the Laplacian
%     lapu = -[(lapux*U(:,1)+lapuy*U(:,2)+lapuz*U(:,3)) (lapvx*U(:,1)+lapvy*U(:,2)+lapvz*U(:,3)) (lapwx*U(:,1)+lapwy*U(:,2)+lapwz*U(:,3))];
% 
%     % Compute the covariant derivative
%     covU = -[(covdux*U(:,1)+covduy*U(:,2)+covduz*U(:,3)) (covdvx*U(:,1)+covdvy*U(:,2)+covdvz*U(:,3)) (covdwx*U(:,1)+covdwy*U(:,2)+covdwz*U(:,3))];
%     covrep = (U(:,1).*(Lx*U(:,1)) + U(:,2).*(Ly*U(:,2)) + U(:,3).*(Lz*U(:,3)));
%     covU = covU + repmat(covrep,1,3);
%     
%     % Compute the coriolis force
%     coriolis = [(-Z*U(:,2) + Y*U(:,3)) (Z*U(:,1) - X*U(:,3)) (-Y*U(:,1) + X*U(:,2))];
%     coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
%     
%     fU = makeGaneshForcing1(N0, x, t, nu, Lx, Ly, Lz, covdux, covduy, covduz, covdvx, covdvy, covdvz, covdwx, covdwy, covdwz);
%     
%     U = -covU + nu*lapu - coriolis + fU;
% U = -covU + nu*lapu;
%     %-----Leray Projection
%     % TODO: build gamdel in a vectorized way.
% 
%     for i=1:N
%         buff(1,:) = d(x(i,:))';
%         buff(2,:) = e(x(i,:))';
%         gamdel(2*i-1:2*i) = buff*U(i,:)';
%     end
%     
%     albet = Adiv\gamdel;
%    
%     % recover the interpolation coefficients
%     interpcoeffs = zeros(N,3);
%     for i = 1:N
%         interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
%     end
%     
%     % TODO:  improve performance
%     Udivfree = zeros(N,3);
%     for i = 1:N
%         for j = 1:N
%             Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
%         end
%     end
%     
%     U = Udivfree;
%     k2 = U;
%     
%     U = Ubuff + (h/2)*k2;
%     %=============================RK4 part 3===============================
%     % Compute the Laplacian
%     lapu = -[(lapux*U(:,1)+lapuy*U(:,2)+lapuz*U(:,3)) (lapvx*U(:,1)+lapvy*U(:,2)+lapvz*U(:,3)) (lapwx*U(:,1)+lapwy*U(:,2)+lapwz*U(:,3))];
% 
%     % Compute the covariant derivative
%     covU = -[(covdux*U(:,1)+covduy*U(:,2)+covduz*U(:,3)) (covdvx*U(:,1)+covdvy*U(:,2)+covdvz*U(:,3)) (covdwx*U(:,1)+covdwy*U(:,2)+covdwz*U(:,3))];
%     covrep = (U(:,1).*(Lx*U(:,1)) + U(:,2).*(Ly*U(:,2)) + U(:,3).*(Lz*U(:,3)));
%     covU = covU + repmat(covrep,1,3);
%     
%     % Compute the coriolis force
%     coriolis = [(-Z*U(:,2) + Y*U(:,3)) (Z*U(:,1) - X*U(:,3)) (-Y*U(:,1) + X*U(:,2))];
%     coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
%     
%     fU = makeGaneshForcing1(N0, x, t, nu, Lx, Ly, Lz, covdux, covduy, covduz, covdvx, covdvy, covdvz, covdwx, covdwy, covdwz);
%     
%     U = -covU + nu*lapu - coriolis + fU;
% U = -covU + nu*lapu;
%     %-----Leray Projection
%     % TODO: build gamdel in a vectorized way.
% 
%     for i=1:N
%         buff(1,:) = d(x(i,:))';
%         buff(2,:) = e(x(i,:))';
%         gamdel(2*i-1:2*i) = buff*U(i,:)';
%     end
%     
%     albet = Adiv\gamdel;
%    
%     % recover the interpolation coefficients
%     interpcoeffs = zeros(N,3);
%     for i = 1:N
%         interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
%     end
%     
%     % TODO:  improve performance
%     Udivfree = zeros(N,3);
%     for i = 1:N
%         for j = 1:N
%             Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
%         end
%     end
%     
%     U = Udivfree;
%     k3 = U;
%     
%     U = Ubuff + h*k3;
%     %=============================RK4 part 4===============================
%     % Compute the Laplacian
%     lapu = -[(lapux*U(:,1)+lapuy*U(:,2)+lapuz*U(:,3)) (lapvx*U(:,1)+lapvy*U(:,2)+lapvz*U(:,3)) (lapwx*U(:,1)+lapwy*U(:,2)+lapwz*U(:,3))];
% 
%     % Compute the covariant derivative
%     covU = -[(covdux*U(:,1)+covduy*U(:,2)+covduz*U(:,3)) (covdvx*U(:,1)+covdvy*U(:,2)+covdvz*U(:,3)) (covdwx*U(:,1)+covdwy*U(:,2)+covdwz*U(:,3))];
%     covrep = (U(:,1).*(Lx*U(:,1)) + U(:,2).*(Ly*U(:,2)) + U(:,3).*(Lz*U(:,3)));
%     covU = covU + repmat(covrep,1,3);
%     
%     % Compute the coriolis force
%     coriolis = [(-Z*U(:,2) + Y*U(:,3)) (Z*U(:,1) - X*U(:,3)) (-Y*U(:,1) + X*U(:,2))];
%     coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
%     
%     fU = makeGaneshForcing1(N0, x, t, nu, Lx, Ly, Lz, covdux, covduy, covduz, covdvx, covdvy, covdvz, covdwx, covdwy, covdwz);
%     
%     U = -covU + nu*lapu - coriolis + fU;
% U = -covU + nu*lapu;
%     %-----Leray Projection
%     % TODO: build gamdel in a vectorized way.
% 
%     for i=1:N
%         buff(1,:) = d(x(i,:))';
%         buff(2,:) = e(x(i,:))';
%         gamdel(2*i-1:2*i) = buff*U(i,:)';
%     end
%     
%     albet = Adiv\gamdel;
%    
%     % recover the interpolation coefficients
%     interpcoeffs = zeros(N,3);
%     for i = 1:N
%         interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
%     end
%     
%     % TODO:  improve performance
%     Udivfree = zeros(N,3);
%     for i = 1:N
%         for j = 1:N
%             Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
%         end
%     end
%     
%     U = Udivfree;
%     k4 = U;
%     
%     %======================Stitch it all together==========================
%     U = Ubuff + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
%     
    %==============================Visualize===============================

    
    %==========================Determine error=============================
    
    
%     Uganesh = makeGaneshTest1(N0, X, t, nu);
%     
%     error = (Uganesh-U).^2;
%     error = sum(sqrt(sum(error,2)));
%     disp(error);
%     U = sqrt(sum((U - Uganesh).^2,2));
end

    function Udiv = getDivFree(U, dmat, emat, Afull, PSIdiv)
        % Get coefficients
        for i = 1:size(U,1)
            gamdel((2*i-1):(2*i)) = [dmat(i,:);emat(i,:)]*(U(i,:).');
        end

        albet = Afull\gamdel;

        coeffs = repmat(albet(1:2:size(albet,1)),1,3).*dmat + ...
                 repmat(albet(2:2:size(albet,1)),1,3).*emat;

        % Finally, recover the div-free part.
        Udiv = PSIdiv*reshape(coeffs',size(coeffs,1)*size(coeffs,2),[]);
        Udiv = reshape(Udiv,3,[])';
        
    end

end