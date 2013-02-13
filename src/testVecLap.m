function [ Maxerr, L2err, kappa ] = testVecLap( x, W, U0, epsilon )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Sense the size of the problem from the initial conditions.  The
    % scalars in U0 are assumed to identify with the R3 vectors in X.  If
    % these are not the same length, function should return error.
    % TODO:  return error
    N = size(U0,1);
    
    % Parameter to all the RBF calls.  This parameter will have to be
    % changed if a different RBF kernel is used.
    twoeps2 = 2*epsilon*epsilon;
   

    % Build the distance matrix.  Simplification is due to the centers
    % being on the surface of the sphere.
    r2 = 2*(1 - x(:,1)*x(:,1).' - x(:,2)*x(:,2).' - x(:,3)*x(:,3).');
    r2 = epsilon*epsilon*r2;
    A = exp(-r2);
    
    kappa = cond(A);


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
 
% attempt 1  
    lapux = -Z*Lz + Y*Y*Lzz - Y*Ly - 2*Y*Z*Lyz + Z*Z*Lyy;
    lapuy = -X*Y*Lzz + X*Z*Lyz + Y*Lx + Y*Z*Lxz - Z*Z*Lxy;
    lapuz =  X*Y*Lyz - X*Z*Lyy + Z*Lx - Y*Y*Lxz + Y*Z*Lxy;

    lapvx = X*Ly - X*Y*Lzz + X*Z*Lyz + Y*Z*Lxz - Z*Z*Lxy;
    lapvy = -X*Lx + X*X*Lzz - Z*Lz - 2*X*Z*Lxz + Z*Z*Lxx;
    lapvz = -X*X*Lyz + X*Y*Lxz + Z*Ly + X*Z*Lxy - Y*Z*Lxx;
    

    lapwx =  X*Lz + X*Y*Lyz - X*Z*Lyy - Y*Y*Lxz + Y*Z*Lxy;
    lapwy =  Y*Lz - X*X*Lyz + X*Y*Lxz + X*Z*Lxy - Y*Z*Lxx;
    lapwz = -Y*Ly + X*X*Lyy - X*Lx - 2*X*Y*Lxy + Y*Y*Lxx;

% attempt 2
%    lapux = -Y*Ly - Y*Y*Lzz + Z*Lz + 2*Y*Z*Lyz - Z*Z*Lyy;
%    lapuy = -Y*Lx + X*Y*Lzz - Z*X*Lyz - Y*Z*Lxz + Z*Z*Lxy;
%    lapuz = -X*Y*Lyz + Y*Y*Lxz - Z*Lx + Z*X*Lyy - Z*Y*Lxy;
    
%    lapvx = -X*Ly + X*Y*Lzz - Z*X*Lyz - Z*Y*Lxz + Z*Z*Lxy;
%    lapvy = X*Lx - X*X*Lzz + Z*Lz + 2*X*Z*Lxz - Z*Z*Lxx;
%    lapvz = X*X*Lyz - X*Y*Lxz - Z*Ly - X*Z*Lxy + Y*Z*Lxx;
    
%    lapwx = -X*Lz - X*Y*Lyz + Y*Y*Lxz + X*Z*Lyy - Y*Z*Lxy;
%    lapwy = X*X*Lyz - Y*Lz - X*Y*Lxz - X*Z*Lxy + Y*Z*Lxx;
%    lapwz = X*Lx - X*X*Lyy + Y*Ly + 2*X*Y*Lxy - Y*Y*Lxx;

% attempt 3
%    lapux = -Y*Z*Lxz + X*Y*Lz + Y*X*X*Lxz + X*Y*Y*Lyz + Z*Lxy - X*Z*Ly - X*X*Z*Lxy - X*Y*Z*Lyy + X*Y*Z*Lzz - X*Z*Z*Lyz;
%    lapuy = Y*Lyz - Y*Y*Lz - X*Y*Y*Lxz - Y*Y*Y*Lyz - Z*Lyy + X*Z*Lx + 2*Y*Z*Ly + X*Y*Z*Lxy + Y*Y*Z*Lyy - Z*Y*Y*Lzz + Z*Z*Lz + Z*Z*Lyz;
%    lapuz = Y*Lzz - X*Y*Lx - Y*Y*Ly - Z*Lyz - 2*Y*Z*Lz - X*Y*Z*Lxz - Y*Y*Z*Lyz + Z*Z*Ly + X*Z*Z*Lxy + Y*Z*Z*Lyy - Y*Z*Z*Lzz - Z*Z*Z*Lyz;

    % Now apply the vector Laplacian to the input data.
    % 
    U = -[(lapux*U0(:,1)+lapuy*U0(:,2)+lapuz*U0(:,3)) (lapvx*U0(:,1)+lapvy*U0(:,2)+lapvz*U0(:,3)) (lapwx*U0(:,1)+lapwy*U0(:,2)+lapwz*U0(:,3))];
    
%    U = [(lapux*U0(:,1)+lapvx*U0(:,2)+lapwx*U0(:,3)) (lapuy*U0(:,1)+lapvy*U0(:,2)+lapwy*U0(:,3)) (lapuz*U0(:,1)+lapvz*U0(:,2)+lapwz*U0(:,3))];

    % Figure out the L-infty residual
    avg = mean(mean(U0./U));
    err = abs(U0 - avg*U);
    errsqrt = sqrt(err(:,1).^2 + err(:,2).^2 + err(:,3).^2);
    Maxerr = max(errsqrt);
    
    % Figure out the L-2 residual
    L2err = 0;
    for i = 1:N
        L2err = L2err + W(i)*(err(i,:)*err(i,:)');
    end
    L2err = sqrt(L2err);
 


end

