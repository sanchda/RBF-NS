function U = rbf_ns(U0, x, h, M, alpha, eps)
% Solves the NS equations on the surface of the sphere using an RBF
% embedded, narrow-band method.  In this method, the differential operators
% are expressed in R3, but restricted to acting on the surface of the
% sphere.  The centers x are assumed to be on the surface of the sphere
% (in a "narrow band" about this surface), thus preserving the action of
% the operators.  These operators (in this case, only the scalar Laplacian)
% have been expressed in terms of differential operators on the surface of
% the sphere (TODO: link to Mathematica notebook or derived .pdf).  The
% differential operators are implemented as RBF differentiation matrices
% following the method described in Fasshauer "RBF Collocation Methods and
% Pseudospectral Methods" (http://www.math.iit.edu/~fass/RBFPS.pdf as of
% 9/22/10).  The kernel used is the Gaussian.

% TODO: lit citation

    % Sense the size of the problem from the initial conditions.  The
    % scalars in U0 are assumed to identify with the R3 vectors in X.  If
    % these are not the same length, function should return error.
    % TODO:  return error
    N = size(U0,1);
    
    % Parameter to all the RBF calls.  This parameter will have to be
    % changed if a different RBF kernel is used.
    twoeps2 = 2*eps*eps;
   

    % Build the distance matrix.  Simplification is due to the centers
   
    % being on the surface of the sphere.
    
    r2 = 2*(1 - x(:,1)*x(:,1).' - x(:,2)*x(:,2).' - x(:,3)*x(:,3).');
    r2 = eps*eps*r2;
    A = exp(-r2);
    
    % Build the matrix for matrix-valued RBF interpolation.
    % This is a gigantic shitmess.  It's not yet obvious how it should be
    % cleaned up.
    e2  = @(x,y) twoeps2*exp((-1 + dot(x,y))*twoeps2);  % Only for S2
    PSI1 = @(x,y) [(-x(2)*y(2) - x(3)*y(3) + twoeps2*(x(3)*y(2) - x(2)*y(3))^2);y(1)*(x(2) + twoeps2*(x(1) - x(3) - y(1) + y(3))*(x(3)*y(2) - x(2)*y(3)));(x(3)*y(1) + twoeps2*(x(2)*y(1) - x(1)*y(2))*(x(3)*y(2) - x(2)*y(3)))]; 
    PSI2 = @(x,y) [(x(1)*(y(2) + twoeps2*(x(1) - x(3) - y(1) + y(3))*(x(3)*y(2) - x(2)*y(3))));(2*x(1)*y(1)*(-1 - (-2 + x(2)*x(2) - 2*x(3)*y(1) + y(2)*y(2) + 2*x(1)*(x(3) + y(1) - y(3)) + twoeps2*(x(3) + y(1))*y(3))));(-x(1)*(-y(2) + twoeps2*(x(2)*y(1) - x(1)*y(2))*(-x(1) + x(3) + y(1) - y(3))))];
    PSI3 = @(x,y) [((x(1)*y(3) - twoeps2*(x(2)*y(1) - x(1)*y(2))*(-x(3)*y(2) + x(2)*y(3))));(-y(1)*(-x(2)+twoeps2*(x(2)*y(1) - x(1)*y(2))*(-x(1) + x(3) + y(1) - y(3))));(-x(1)*y(1)-x(2)*y(2)+twoeps2*(x(2)*y(1)-x(1)*y(2))^2)];
    PSI  = @(x,y) e2(x,y)*[PSI1(x,y)';PSI2(x,y)';PSI3(x,y)'];
    
    dx   = @(x) (1/sqrt(1 - x(3)^2))*[(-x(3)*x(1)) (x(3)*x(2)) (1-x(3)*x(3))]';
    ex   = @(x) (1/sqrt(1 - x(3)^2))*[-x(2) x(1) 0]';
    
    
    mat1 = @(i,j) [-dot(ex(x(i,:)),ex(x(j,:))) dot(ex(x(i,:)),dx(x(j,:))); dot(dx(x(i,:)),ex(x(j,:))) -dot(dx(x(i,:)),dx(x(j,:)))];
    mat2 = @(i,j) [dot(ex(x(i,:)),x(j,:)) + dot(ex(x(j,:)),x(i,:))  dot(ex(x(i,:)),x(j,:)) - dot(dx(x(j,:)),x(i,:));-dot(dx(x(i,:)),x(j,:)) + dot(ex(x(j,:)),x(i,:)) -dot(dx(x(i,:)),x(j,:)) - dot(dx(x(j,:)),x(i,:))];
    eta1 = @(r) -twoeps2*exp(-r*r*eps*eps);
    mu1  = @(r) twoeps2*twoeps2*exp(-r*r*eps*eps);
    r    = @(i,j) sqrt( (x(i,1) - x(j,1))^2 + (x(i,2) - x(j,2))^2 + (x(i,3) - x(j,3))^2);
    eta  = @(i,j) eta1(r(i,j));
    mu   = @(i,j) mu1(r(i,j));
    
    A2   = @(i,j) eta(i,j)*mat1(i,j) + mu(i,j)*mat2(i,j);
    
    % TODO: build Adiv in a more intelligent way
    Adiv = zeros(2*N,2*N);
    
    for i = 1:N
        for j = 1:N
            Adiv((2*i-1):(2*i),(2*j-1):(2*j)) = A2(i,j);
        end
    end
    
    % TODO: build gamdel in a better way.  Also, seriously--did you wake
    % up on the retarded side of the bed when you named this array?
    gamdel = zeros(2*N,1);
    buff = zeros(2,3);

    for i=1:N
        buff(1,:) = dx(x(i,:))';
        buff(2,:) = ex(x(i,:))';
        gamdel(2*i-1:2*i) = buff*U0(i,:)';
    end
    
    albet = ((gamdel')/Adiv)';
    
    % Now, recover the coefficients.
    PSImat = zeros(3*N,3*N);
    PSIarg = zeros(3*N,1);
    
    for i = 1:N
        for j = 1:N
            PSImat((3*i-2):(3*i),(3*j-2):(3*j)) = PSI(x(i,:),x(j,:));
        end
        PSIarg((3*i-2):(3*i),:) = albet(2*i-1)*dx(x(i,:)) + albet(2*i)*ex(x(i,:));
    end
    
    Udivbuff = PSImat*PSIarg;
    
    Udivfree = reshape(Udivbuff,N,3);
    
    U = Udivfree;
    
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
    

    % Initialize the vector Laplacian.
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
  
    lapux = -Z*Lz + Y*Y*Lzz - Y*Ly - 2*Y*Z*Lyz + 2*Z*Z*Lyy;
    lapuy = -X*Y*Lzz + X*Z*Lyz + Y*Lx + Y*Z*Lxz - Z*Z*Lxy;
    lapuz =  X*Y*Lyz - X*Z*Lyy + Z*Lx - Y*Y*Lxz + Y*Z*Lxy;
    
    lapvx =  X*Ly + X*Z*Lyz - X*Y*Lyy + Y*Z*Lxz + Z*Z*Lxy;
    lapvy = -Z*Lz + X*X*Lzz - X*Lx - 2*X*Z*Lx + Z*Z*Lxx;
    lapvz =  Z*Ly - X*X*Lyz + X*Y*Lxz + X*Z*Lxy - Y*Z*Lxx;
    
    lapwx =  X*Lz + X*Y*Lyz - X*Z*Lyy + Y*Y*Lxz + Y*Z*Lxy;
    lapwy =  Y*Lz + X*X*Lyz + X*Y*Lxz + X*Z*Lxy - Y*Z*Lxx;
    lapwz = -Y*Ly + X*X*Lyy - X*Lx + 2*X*Y*Lxy + Y*Y*Lxx;
       
    % TODO: multiply the result of lap(U) by mu.
    
    
    % Clear up the unused matrices
%    clear('X','Y','Z');
%    clear('Lx','Ly','Lz');
%    clear('Lxy','Lxz','Lyz');
%    clear('Lxx','Lyy','Lzz');
 
	% Get the first timestep
    U = U0;
 
	% Now commence timestepping.  This is the standard RK4 method, which
	% may or may not be well-suited for this problem.
    for i=1:M
        disp(i);
        rk1 = h*rkfun(U, lapu);
        rk2 = h*rkfun(U + (1/2)*rk1, lapu);
        rk3 = h*rkfun(U + (1/2)*rk2, lapu);
        rk4 = h*rkfun(U + rk3, lapu);
        
        U = U + (1/6)*(rk1 + 2*rk2 + 2*rk3 + rk4);
    end
    
    function Un = rkfun(U, lapu)
       Un = alpha*lapu*U;
    end

    % TODO: visualization, analysis, etc.  Currently on the final timestep
    % is visualized from the driver.

end