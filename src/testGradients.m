%% Initialize!
cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';

%==========================================================================
%                         Parameters and Constants                        
%==========================================================================
N = 31;        % Somehow related to the number of centers.  For the
               % ME points, the number of centers is (N+1)^2.

% Generated parameters
divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);
eps_Leray = divFree_geteps(N);
beta = 12;
c = (4*pi)^2;
eps_PDE = (beta/c)*(N+1)^(9/8);

disp('Set parameters')
%==========================================================================
%                            Generate RBF nodes                        
%==========================================================================
x = getMEPoints(N);
W = x(:,4);
x = x(:,1:3);
% Rotate X through by a small angle to avoid poles
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    x(i,:) = (theta*x(i,:)')';
end
x = sortrows(x,3);

disp('RBF nodes loaded')


%==========================================================================
%                         Setup Diff. Operators
%==========================================================================  
twoeps2 = 2*eps_PDE*eps_PDE;

% Scalar RBF kernel and its derivative
syms r;

phi = @(r2) exp(-eps_PDE*eps_PDE*r2);

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

%====================Initialize Projection Operator=============
% explicit indexing in case numel(x)>3 for some reason
P = @(x) eye(3) - [x(:,1);x(:,2);x(:,3)]*[x(:,1) x(:,2) x(:,3)];


% Make a sparse, block-diagonal matrix with copies of P.
Pxmat = zeros(3*size(x,1),3*size(x,1));
for i = 1:size(x,1)
   %TODO:  built in a sparse way?
   Pxmat((3*i-2):(3*i),(3*i-2):(3*i)) = P(x(i,:));
end
Pxmat = sparse(Pxmat);
    


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

grad1 = [gradx;grady;gradz];
grad2 = [Lx; Ly; Lz];
grad2 = Pxmat*grad2;

disp('Differential operators created')

%% Test Gradient Operators

U = sin(x(:,1)).*cos(x(:,1)) + x(:,1).*x(:,2).^2 + x(:,1).*exp(x(:,1));
    
grad1U = reshape(grad1*U,[],3);
grad2U = reshape(grad2*U,[],3);



gradUR3 = reshape(gradR3*U,[],3);
gradUR3 = reshape(Pxmat*reshape(gradUR3',[],1),3,[])';

err = grad1U - grad2U;

err = sqrt(err(:,1).^2 + err(:,2).^2 + err(:,3).^2);
plot(err)

%% Test Vector Laplacian