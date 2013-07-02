% This is a test to reproduce the following bug:
% Pxmat should project onto the tangent space on the sphere, but doesn't in
% at least this case.
%
% Operating hypotheses:
% Pxmat is idempotent(tested; see case 4)
% Pxmat projects onto tangent space (hypothesis)

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
U  = getDivFree(2,x); 
U  = U(:,1:3);

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

%% Case 2:  projection check:  projection of normals is tiny
% Note that the resulting field is so small the normalized visualization is
% pointless
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
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 3a:  behavior of algebraic synthesis with Px and grad (i.e., gradS2)
% The result does not appear restricted to the tangent space
U = x;
 
% Just pick a direction, since grad accepts only vectors.
U = gradS2*U(:,1);
U = reshape(U,[],3);
 
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
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
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off

%% Case 4:  pxmat is idempotent
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
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off
