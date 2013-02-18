function U = navierstokes(x, U0, H, h, M, epsilon, nu, omega)
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
    
    
%==========================================================================
%                             Setup Leray Projector
%==========================================================================

    % Zonal and meridional bases
    d = @(x) (1/sqrt(1-x(3)^2))*[(-x(3)*x(1)) (-x(3)*x(2)) (1-x(3)^2)];
    e = @(x) (1/sqrt(1-x(3)^2))*[(-x(2)) x(1) 0];
    
    % The function Q(x) generates the matrix which will project a vector in R3
    % onto the tangent space of the sphere at x \in S2
    Q = @(x) [0 x(3) (-x(2)); (-x(3)) 0 x(1);  x(2) (-x(1)) 0];
    
    % P(x) projects into the normal space
    P = @(x) eye(3) - [x(1);x(2);x(3)]*[x(1) x(2) x(3)];
    
    % Construct the matrix RBF
    PSIdiv  = @(x,y) (Q(x)')*(-H(x,y,epsilon))*Q(y);
%    PSIcrl  = @(x,y) (P(x)')*(-H(x,y,epsilon))*P(y);
    PSI     = @(x,y) PSIdiv(x,y) + PSIcrl(x,y);
    
    Adiv =  makeAdiv(x, PSI, d, e);
    
%==========================================================================
%                         Setup Diff. Operators
%==========================================================================  
    
    % Parameter to all the RBF calls.  This parameter will have to be
    % changed if a different RBF kernel is used.
    twoeps2 = 2*epsilon*epsilon;

    % Build the distance matrix.  Simplification is due to the centers
    % being on the surface of the sphere.
    r2 = 2*(1 - x(:,1)*x(:,1).' - x(:,2)*x(:,2).' - x(:,3)*x(:,3).');
    r2 = epsilon*epsilon*r2;
    A = exp(-r2);

    % Initialize the differentiation matrices
    % ASSERT:  X contains only points on the surface of the unit sphere
    
    % In order to write things like (x1-y1)*A (see accompanying Mathematica
    % notebook), we need to form some distance matrices.  The following
    % method was inspired by Joseph Kirk, author of distmat.
    % (http://www.mathworks.com/matlabcentral/fileexchange/15145-distance-m
    % atrix/content/distmat.m case 2)

    % TODO: time whether a similar method would yield more efficient
    % initialization for A.
    xdiff1 = reshape(x(:,1),1,N,1);
    xdiff2 = reshape(x(:,1),N,1,1);
    
    ydiff1 = reshape(x(:,2),1,N,1);
    ydiff2 = reshape(x(:,2),N,1,1);
    
    zdiff1 = reshape(x(:,3),1,N,1);
    zdiff2 = reshape(x(:,3),N,1,1);
    
    xdiff = -xdiff1(ones(N,1),:,:) + xdiff2(:,ones(N,1),:);
    ydiff = -ydiff1(ones(N,1),:,:) + ydiff2(:,ones(N,1),:);
    zdiff = -zdiff1(ones(N,1),:,:) + zdiff2(:,ones(N,1),:);
    
    % To form Lx, a matrix differential operator, the same process as
    % before is used, but where the RBF kernel is the derivative (in x)
    % of the one before.  Since this function is the exponential, we can
    % reuse some of the previous data in the derivative (this intermediate
    % matrix should really be called something like Bx, but we reuse Lx
    % here).  Finally, Lx = Bx*inv(A), but inv() isn't stable so we do a
    % right-division.
   
    Lx = -twoeps2*A.*xdiff;
    Lx = Lx/A;

    Ly = -twoeps2*A.*ydiff;  	
    Ly = Ly/A;
 	
    Lz = -twoeps2*A.*zdiff;  	
    Lz = Lz/A;
    
    
    Lxy = twoeps2*twoeps2*A.*xdiff.*ydiff;
    Lxy = Lxy/A;
  	
    Lxz = twoeps2*twoeps2*A.*xdiff.*zdiff;
    Lxz = Lxz/A;
    
    Lyz = twoeps2*twoeps2*A.*ydiff.*zdiff;
    Lyz = Lyz/A;

    
    Lxx = twoeps2*A.*(-1 + twoeps2*(xdiff.*xdiff)); 	
    Lxx = Lxx/A;	  	
    
    Lyy = twoeps2*A.*(-1 + twoeps2*(ydiff.*ydiff));  	
    Lyy = Lyy/A;

    Lzz = twoeps2*A.*(-1 + twoeps2*(zdiff.*zdiff));
    Lzz = Lzz/A;

    
    % Clear out the unused matrices
    clear('xdiff1','xdiff2','xdiff');
    clear('ydiff1','ydiff2','ydiff');
    clear('zdiff1','zdiff2','zdiff');
    

%=====================Initialize the vector Laplacian=====================

    % ASSERT:  X contains only points on the surface of the unit sphere
    %
    % Since U = (u,v,w) is a vector, we write veclapU in terms of the
    % action of the vector Laplacian on each individual component of U.
    % Moreover, since the action of each component relies on contributions
    % from each coordinate, we write the three componentwise Laplacians
    % in three parts, for a total of nine operators.
    
    X = diag(x(:,1));
    Y = diag(x(:,2));
    Z = diag(x(:,3));
    
    % Need this later for the coriolis force.  Let's just define it now.
    zsqrt = 1 - Z.*Z;
    zsqrt = sqrt(zsqrt);
 
    lapux = -Z*Lz + Y*Y*Lzz - Y*Ly - 2*Y*Z*Lyz + Z*Z*Lyy;
    lapuy = -X*Y*Lzz + X*Z*Lyz + Y*Lx + Y*Z*Lxz - Z*Z*Lxy;
    lapuz =  X*Y*Lyz - X*Z*Lyy + Z*Lx - Y*Y*Lxz + Y*Z*Lxy;

    lapvx = X*Ly - X*Y*Lzz + X*Z*Lyz + Y*Z*Lxz - Z*Z*Lxy;
    lapvy = -X*Lx + X*X*Lzz - Z*Lz - 2*X*Z*Lxz + Z*Z*Lxx;
    lapvz = -X*X*Lyz + X*Y*Lxz + Z*Ly + X*Z*Lxy - Y*Z*Lxx;  

    lapwx =  X*Lz + X*Y*Lyz - X*Z*Lyy - Y*Y*Lxz + Y*Z*Lxy;
    lapwy =  Y*Lz - X*X*Lyz + X*Y*Lxz + X*Z*Lxy - Y*Z*Lxx;
    lapwz = -Y*Ly + X*X*Lyy - X*Lx - 2*X*Y*Lxy + Y*Y*Lxx;
    
    
%=====================Initialize covariant Derivative======================

    % Reuses information from the above section 
    covdux = X*Y*Lz - X*Z*Ly;
    covduy = -X*X*Lz + X*Z*Lx;
    covduz = X*X*Ly - X*Y*Lx;
    
    covdvx = Y*Y*Lz - Y*Z*Ly;
    covdvy = -X*Y*Lz + Y*Z*Lx;
    covdvz = X*Y*Ly - Y*Y*Lx;
    
    covdwx = Y*Z*Lz - Z*Z*Ly;
    covdwy = -X*Z*Lz + Z*Z*Lx;
    covdwz = X*Z*Ly - Y*Z*Lx;
    
    
%==========================================================================
%                             Begin Timestepping
%==========================================================================  
 
% pre-defines
U = U0;
gamdel = zeros(2*N,1);
buff = zeros(2,3);
dmat = zeros(N,3);
emat = zeros(N,3);

for i = 1:N
    dmat(i,:) = d(x(:,1));
    emat(i,:) = e(x(:,1));
end

for c = 1:M

    %=======================Solve approximate system=======================
    
    % Compute the Laplacian
    lapu = -[(lapux*U(:,1)+lapuy*U(:,2)+lapuz*U(:,3)) (lapvx*U(:,1)+lapvy*U(:,2)+lapvz*U(:,3)) (lapwx*U(:,1)+lapwy*U(:,2)+lapwz*U(:,3))];

    % Compute the covariant derivative
    covU = -[(covdux*U(:,1)+covduy*U(:,2)+covduz*U(:,3)) (covdvx*U(:,1)+covdvy*U(:,2)+lapvz*U(:,3)) (covdwx*U(:,1)+covdwy*U(:,2)+covd*U(:,3))];
    covrep = (U(:,1).*Lx*(U(:,1)) + U(:,2).*Ly*(U(:,2)) + U(:,3).*Ly*(U(:,3)));
    covU = covU + repmat(covrep,1,3);
    
    % Compute the coriolis force
    coriolis = [(-Z*U(:,2) + Y*U(:,3)) (Z*U(:,1) - X*U(:,3)) (-Y*U(:,1) + X*U(2,:))];
    coriolis = 2*omega*repmat(zsqrt,1,3).*coriolis;
    
    fU = -covU + nu*lapu - coriolis + f;
    U = U + h*fU;
       
    %========================Project onto div-free=========================
    
    % TODO: build gamdel in a vectorized way.

    for i=1:N
        buff(1,:) = d(x(i,:))';
        buff(2,:) = e(x(i,:))';
        gamdel(2*i-1:2*i) = buff*U(i,:)';
    end
    
    albet = Adiv\gamdel;
   
    % recover the interpolation coefficients
    interpcoeffs = zeros(N,3);
    for i = 1:N
        interpcoeffs(i,:) = albet(2*i-1)*d(x(i,:)) + albet(2*i)*e(x(i,:));
    end
    
    % TODO:  improve performance
    Udivfree = zeros(N,3);
    for i = 1:N
        for j = 1:N
            Udivfree(i,:) = Udivfree(i,:) + (PSIdiv(x(i,:),x(j,:))*(interpcoeffs(j,:)'))';
        end
    end
    
    %==========================Determine error=============================
   
end
    


end