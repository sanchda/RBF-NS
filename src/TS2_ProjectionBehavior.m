% Test suite for checking various features of the code.  This test suite is
% currently intended for developers and should not be exposed to users.

% TODO:
% Develop quantitative tests for Laplacian, gradient, and covariant
% derivative.
% Use the initialization and building routines utilized by navierstokes.m
% itself--until then WARNING:  routines used herein may not reflect those
% used in live code!

cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';

%% Get cell centers
N = 12;
epsPDE  = 3;
surfeps = 3;
epsilon = epsPDE; % legacy define
twoeps2 = 2*epsPDE*epsPDE;


x = getMEPoints(N);
W = x(:,4); % Keep the quadrature weights
x = x(:,1:3);
% Rotate X through by a small angle to avoid poles
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    x(i,:) = (theta*x(i,:)')';
end

% Debug line--sort X by Z position
x = sortrows(x,3);

%% Generate projection matrix

% explicit indexing in case numel(x)>3 for some reason
P = @(x) eye(3) - [x(:,1);x(:,2);x(:,3)]*[x(:,1) x(:,2) x(:,3)];


% Make a sparse, block-diagonal matrix with copies of P.
Pxmat = zeros(3*size(x,1),3*size(x,1));
for i = 1:size(x,1)
   %TODO:  built in a sparse way?
   Pxmat((3*i-2):(3*i),(3*i-2):(3*i)) = P(x(i,:));
end
Pxmat = sparse(Pxmat);

%% Compute R3 derivative matrices
% Scalar RBF kernel
% ASSERT:  phi is Gaussian; else L? definitions below will fail since the
% exponent is its own derivative
% TODO:  build derivatives in a kernel-agnostic way
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

%% Define the gradient in R3
gradR3 = [Lx;Ly;Lz];

%% Define the gradient in S2
X = diag(x(:,1));
Y = diag(x(:,2));
Z = diag(x(:,3));

gradx = (eye(size(X,1))-X*X)*Lx - X*Y*Ly - X*Z*Lz;
grady = -X*Y*Lx + (eye(size(X,1))-Y*Y)*Ly - Y*Z*Lz;
gradz = -X*Z*Lx - Y*Z*Ly + (eye(size(X,1))-Z*Z)*Lz;

gradS2 = [gradx;grady;gradz];

%% Define the vector Laplacian in S2

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


%% Initialize visualization parameters
%-q suffix denotes it is for the quiver arrows; otherwise, for surface
%colormap interpolation
eps = 1.0;
epsq = 0.1;
clear @phi;
phi = @(r2) exp(-eps*eps*r2);
phiq = @(r2) exp(-epsq*epsq*r2);

s2 = @(t,p) [cos(t).*sin(p) sin(t).*sin(p) cos(p)];
sz = [101, 201];
szq = [25 25];
MM = prod(sz);  
MMq = prod(szq);  
NN = size(x,1);

% surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi,pi,sz(1)));
[llq,ttq]=meshgrid(linspace(-pi,pi,szq(2)),linspace(-pi,pi,szq(1)));
xx=s2(ll(:),tt(:));
xxq=s2(llq(:),ttq(:));

% Interpolate to the grid
re2=(repmat(xx(:,1),[1 NN])-repmat(x(:,1).',[MM 1])).^2;
re2=re2+(repmat(xx(:,2),[1 NN])-repmat(x(:,2).',[MM 1])).^2;
re2=re2+(repmat(xx(:,3),[1 NN])-repmat(x(:,3).',[MM 1])).^2;

re2q=(repmat(xxq(:,1),[1 NN])-repmat(x(:,1).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,2),[1 NN])-repmat(x(:,2).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,3),[1 NN])-repmat(x(:,3).',[MMq 1])).^2;
yy=reshape(xx(:,2),sz);
zz=reshape(xx(:,3),sz);
xx=reshape(xx(:,1),sz);

yyq=reshape(xxq(:,2),szq);
zzq=reshape(xxq(:,3),szq);
xxq=reshape(xxq(:,1),szq);


%% Case 1a:  visualization check:  normal vectors
% Clear current figure
clf;

% Note that the colormap is normalized.
U = x;
max(max(U))
min(min(U))

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 1b:  visualization check:  vector spherical harmonics are tangent
% Clear current figure
clf;

U  = getDivFree(2,x); 
U  = U(:,1:3);

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 2:  projection check:  projection of normals is tiny
% Note that the resulting field is so small the normalized visualization is
% pointless

% Clear current figure
clf;

U = x;

U = Pxmat*reshape(U',[],1);
U = reshape(U',3,[])';

max(max(U))
min(min(U))

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 3a:  behavior of algebraic synthesis with Px and grad (i.e., gradS2)
% Clear current figure
clf;

% The result does not appear restricted to the tangent space
U = x;
 
% Just pick a direction, since grad accepts only vectors.
U = gradS2*U(:,1);
U = reshape(U,[],3);
 
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz);
uu=reshape(phi(re2)*(Achol\(Achol.'\U(:,1))),sz);
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 
 
% Plot the results
hold on
% quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 3b:  behavior of projected gradR3
% Clear current figure
clf;

U = x;

% Just pick a direction, since grad accepts only vectors.
U = gradR3*U(:,1);
U = reshape(U,[],3);
U = Pxmat*reshape(U',[],1);
U = reshape(U',3,[])';

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 3c:  behavior of lap
% Clear current figure
clf;

U = x;

U = 

% Just pick a direction, since grad accepts only vectors.
lapU = lap*reshape(U,[],1);
lapU = reshape(lapU,[],3);
U = lapU;

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off


%% Case 4:  pxmat is idempotent
% Clear current figure
clf;

U = x;

% Just pick a direction, since grad accepts only vectors.
U = gradR3*U(:,1);
U = reshape(U,[],3);
U = Pxmat*reshape(U',[],1);
U = reshape(U',3,[])';
U = reshape(U,[],3);
U = Pxmat*reshape(U',[],1);
U = reshape(U',3,[])';

uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(x(:,1),x(:,2),x(:,3),U(:,1),U(:,2),U(:,3),1)
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off
