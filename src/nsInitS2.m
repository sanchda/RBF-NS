function [lap grad Lx Ly Lz Achol Afull Acrl Adiv PSIfullmat PSIcrlmat PSIdivmat Pxmat] = nsInitS2(x, H, epsilon)
% Initializes everything needed by the Navier-Stokes spherical code.

%==========================================================================
%                               Initialization
%==========================================================================
    % Sense the size of the problem from the initial conditions.  The
    % scalars in U0 are assumed to identify with the R3 vectors in X.  If
    % these are not the same length, function should return error.
    % TODO:  return error
    N = size(x,1);
    twoeps2 = 2*epsilon*epsilon;
    
    
%==========================================================================
%                             Setup Leray Projector
%==========================================================================

    % Efficient element-wise multiplication of the columns of A against the
    % vector v.
    timesMatVec = @(A,b) bsxfun(@times,A,b(:));

    % Zonal and meridional bases.  Works on column-vectors as well as
    % arrays.
    d = @(x) timesMatVec([(-x(:,3).*x(:,1)) (-x(:,3).*x(:,2)) ...
        (1-x(:,3).^2)],(1./sqrt(1-x(:,3).^2)));
    
    e = @(x) timesMatVec([(-x(:,2)) x(:,1) ...
        0*x(:,1)],(1./sqrt(1-x(:,3).^2)));
    
    % The function Q(x) generates the matrix which will project a vector
    % in R3 onto the tangent space of the sphere at x in S2
    Q = @(x) [0*x(:,1) x(:,3) (-x(:,2));...
              (-x(:,3)) 0*x(:,1) x(:,1);...
                x(:,2) (-x(:,1)) 0*x(:,3)];
            
    % explicit indexing in case numel(x)>3 for some reason
    P = @(x) eye(3) - [x(:,1);x(:,2);x(:,3)]*[x(:,1) x(:,2) x(:,3)];
    
    
    % Make a sparse, block-diagonal matrix with copies of P.
    Pxmat = zeros(3*size(x,1),3*size(x,1));
    for i = 1:size(x,1)
       %TODO:  built in a sparse way?
       Pxmat((3*i-2):(3*i),(3*i-2):(3*i)) = P(x(i,:));
    end
    Pxmat = sparse(Pxmat);
    
    % Construct the matrix RBF
    PSIdiv  = @(x,y) (Q(x)')*(-H(x,y,epsilon))*Q(y);
    PSIcrl  = @(x,y) (P(x)')*(-H(x,y,epsilon))*P(y);
    PSI     = @(x,y) PSIdiv(x,y) + PSIcrl(x,y);
    
    Afull =  makeSBFKernel(x, PSI, d, e);
    Acrl =   makeSBFKernel(x, PSIcrl, d, e);
    Adiv =   makeSBFKernel(x, PSIdiv, d, e);
    
    PSIfullmat = zeros(3*size(x,1),3*size(x,1));
    PSIdivmat = PSIfullmat;
    PSIcrlmat = PSIfullmat;
    
    for i=1:size(x,1)
        for j=1:size(x,1)
            PSIfullmat((3*i-2):(3*i),(3*j-2):(3*j)) = PSI(x(i,:),x(j,:));
            PSIdivmat((3*i-2):(3*i),(3*j-2):(3*j)) = PSIdiv(x(i,:),x(j,:));
            PSIcrlmat((3*i-2):(3*i),(3*j-2):(3*j)) = PSIcrl(x(i,:),x(j,:));
        end
    end
    
    
%==========================================================================
%                         Setup Diff. Operators
%==========================================================================  
    
    % Scalar RBF kernel and its derivative
    syms r;

    phi = @(r2) exp(-epsilon*epsilon*r2);

    % Build the distance matrix.
    xdist = repmat(x(:,1),[1 size(x,1)]);
    xdist = xdist - xdist.';
    ydist = repmat(x(:,2),[1 size(x,1)]);
    ydist = ydist - ydist.';
    zdist = repmat(x(:,3),[1 size(x,1)]);
    zdist = zdist - zdist.';
    
    A = phi(xdist.^2 + ydist.^2 + zdist.^2);
    Achol = chol(A);

    % Initialize the differentiation matrices
    % ASSERT:  X contains only points on the surface of the unit sphere
    
    % To form Lx, a matrix differential operator, the same process as
    % before is used, but where the RBF kernel is the derivative (in x)
    % of the one before.  Since this function is the exponential, we can
    % reuse some of the previous data in the derivative (this intermediate
    % matrix should really be called something like Bx, but we reuse Lx
    % here).  Finally, Lx = Bx*inv(A), but inv() isn't stable so we do a
    % right-division.
   
    Lx = -twoeps2*A.*xdist;
    Lx = (Lx/Achol)/Achol.';

    Ly = -twoeps2*A.*ydist;  	
    Ly = (Ly/Achol)/Achol.';
 	
    Lz = -twoeps2*A.*zdist;  	
    Lz = (Lz/Achol)/Achol.';
    
    
    Lxy = twoeps2*twoeps2*A.*xdist.*ydist;
    Lxy = (Lxy/Achol)/Achol.';
  	
    Lxz = twoeps2*twoeps2*A.*xdist.*zdist;
    Lxz = (Lxz/Achol)/Achol.';
    
    Lyz = twoeps2*twoeps2*A.*ydist.*zdist;
    Lyz = (Lyz/Achol)/Achol.';

    
    Lxx = twoeps2*A.*(-1 + twoeps2*(xdist.*xdist)); 	
    Lxx = (Lxx/Achol)/Achol.';	  	
    
    Lyy = twoeps2*A.*(-1 + twoeps2*(ydist.*ydist));  	
    Lyy = (Lyy/Achol)/Achol.';

    Lzz = twoeps2*A.*(-1 + twoeps2*(zdist.*zdist));
    Lzz = (Lzz/Achol)/Achol.';

 
    

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
 
    lapux = -Z*Lz + Y*Y*Lzz - Y*Ly - 2*Y*Z*Lyz + Z*Z*Lyy;
    lapuy = -X*Y*Lzz + X*Z*Lyz + Y*Lx + Y*Z*Lxz - Z*Z*Lxy;
    lapuz =  X*Y*Lyz - X*Z*Lyy + Z*Lx - Y*Y*Lxz + Y*Z*Lxy;

    lapvx = X*Ly - X*Y*Lzz + X*Z*Lyz + Y*Z*Lxz - Z*Z*Lxy;
    lapvy = -X*Lx + X*X*Lzz - Z*Lz - 2*X*Z*Lxz + Z*Z*Lxx;
    lapvz = -X*X*Lyz + X*Y*Lxz + Z*Ly + X*Z*Lxy - Y*Z*Lxx;  

    lapwx =  X*Lz + X*Y*Lyz - X*Z*Lyy - Y*Y*Lxz + Y*Z*Lxy;
    lapwy =  Y*Lz - X*X*Lyz + X*Y*Lxz + X*Z*Lxy - Y*Z*Lxx;
    lapwz = -Y*Ly + X*X*Lyy - X*Lx - 2*X*Y*Lxy + Y*Y*Lxx;
    
    lap = [lapux lapuy lapuz; lapvx lapvy lapvz; lapwx lapwy lapwz];
    
    
    
%=====================Initialize Gradient======================

    % Reuses information from the above section 
    gradx = (eye(size(X,1))-X*X)*Lx - X*Y*Ly - X*Z*Lz;
    grady = -X*Y*Lx + (eye(size(X,1))-Y*Y)*Ly - Y*Z*Lz;
    gradz = -X*Z*Lx - Y*Z*Ly + (eye(size(X,1))-Z*Z)*Lz;
    
    grad = [gradx;grady;gradz];
    
    
%     covux = X*Y*Lz - X*Z*Ly;
%     covuy = -X*X*Lz + X*Z*Lx;
%     covuz = X*X*Ly - X*Y*Lx;
%     
%     covvx = Y*Y*Lz - Y*Z*Ly;
%     covvy = -X*Y*Lz + Y*Z*Lx;
%     covvz = X*Y*Ly - Y*Y*Lx;
%     
%     covwx = Y*Z*Lz - Z*Z*Ly;
%     covwy = -X*Z*Lz + Z*Z*Lx;
%     covwz = X*Z*Ly - Y*Z*Lx;
%     
%     cov = [covux covuy covuz; covvx covvy covvz; covwx covwy covwz];

    
end

