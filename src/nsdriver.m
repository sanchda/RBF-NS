cd 'C:\Users\david\Documents\GitHub\RBF-NS\src';

%==========================================================================
%                         Parameters and Constants                        
%==========================================================================

nu = 1;     % Parameter for the NS equation
omega=0;    % Strength of coriolis force
N = 12;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
N0 = 3;     % Highest spherical harmonic in the test
M = 1;      % how many iterations to run the simulation for
h = 1/(2*N);   % timestep
divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);
epsLeray = divFree_geteps(N);
epsPDE   = 3;
surfeps = 3;

%==============================debug parameters============================
% Check the version of the vector Laplacian saved in testVecLap.m
test_lap = 0;

%==========================================================================
%                            Generate RBF nodes                        
%==========================================================================

X = getMEPoints(N);
W = X(:,4);
X = X(:,1:3);
% Rotate X through by a small angle to avoid poles
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    X(i,:) = (theta*X(i,:)')';
end
X = sortrows(X,3);


%==========================================================================
%                            Hessian for matrix RBF                        
%==========================================================================
HGA =  @(x,y,eps) [exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
                    -1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(1)+(-1).*y(1)).^2),4.*exp( ...
                    1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).* ...
                    y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(2)+(-1).*y(2)),4.*exp(1).^((-1).* ...
                    eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).* ...
                    eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3));4.*exp(1).^((-1).*eps.^2.*(( ...
                    x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(1)+( ...
                    -1).*y(1)).*(x(2)+(-1).*y(2)),exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+( ...
                    x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(2)+( ...
                    -1).*y(2)).^2),4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1) ...
                    .*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)); ...
                    4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
                    -1).*y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3)),4.*exp(1).^(( ...
                    -1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)) ...
                    .*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)),exp(1).^((-1).*eps.^2.*(( ...
                    x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).* ...
                    eps.^2+4.*eps.^4.*(x(3)+(-1).*y(3)).^2)];

[lap grad Lx Ly Lz Achol Afull Acrl Adiv PSIfull PSIcrl PSIdiv Pxmat] = nsInitS2(X, HGA, 1.05*epsLeray);
%==========================================================================
%                     Optional Test: Vector Laplacian                         
%==========================================================================
if test_lap == 1
    
   divFree_L = 1;                                     % how high is the degree of the VSH
   divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);  % empirically fitted shape parameter for the vector Laplacian
   divFree_U0  = getDivFree(2,X); 
   divFree_U0  = divFree_U0(:,1:3);                   % use a divergence-free VSH
   divFree_eps = divFree_geteps(N);
   [~, divFree_Maxerr, divFree_L2err, divFree_kappa] = testVecLap(X, W, divFree_U0, divFree_eps);
   [divFree_L2err divFree_Maxerr]
   
end


%==========================================================================
%                            Generate initial VF                       
%==========================================================================
U0 = makeGaneshTest1(1, X, t, nu);
U0  = getDivFree(2,X); 
U0  = U0(:,1:3);
U = U0;
%==========================================================================
%                        Initialize Visualization Stuff                       
%==========================================================================
% 
% Parametrized surface and constants.  -q suffix denotes it's for the
% quiver plot
eps = 1.0;
epsq = 0.1;
phi = @(r2) exp(-eps*eps*r2);
phiq = @(r2) exp(-epsq*epsq*r2);

s2 = @(t,p) [cos(t).*sin(p) sin(t).*sin(p) cos(p)];
sz = [101, 201];
szq = [25 25];
MM = prod(sz);  
MMq = prod(szq);  
NN = size(X,1);

% surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi,pi,sz(1)));
[llq,ttq]=meshgrid(linspace(-pi,pi,szq(2)),linspace(-pi,pi,szq(1)));
xx=s2(ll(:),tt(:));
xxq=s2(llq(:),ttq(:));

% Interpolate to the grid
re2=(repmat(xx(:,1),[1 NN])-repmat(X(:,1).',[MM 1])).^2;
re2=re2+(repmat(xx(:,2),[1 NN])-repmat(X(:,2).',[MM 1])).^2;
re2=re2+(repmat(xx(:,3),[1 NN])-repmat(X(:,3).',[MM 1])).^2;

re2q=(repmat(xxq(:,1),[1 NN])-repmat(X(:,1).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,2),[1 NN])-repmat(X(:,2).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,3),[1 NN])-repmat(X(:,3).',[MMq 1])).^2;
yy=reshape(xx(:,2),sz);
zz=reshape(xx(:,3),sz);
xx=reshape(xx(:,1),sz);

yyq=reshape(xxq(:,2),szq);
zzq=reshape(xxq(:,3),szq);
xxq=reshape(xxq(:,1),szq);


% Simulate!
for c = 1:1
    
U = navierstokes(X, U0, h, 2, eps, nu, omega, lap, grad, surfeps, Lx, Ly, Lz, Afull, Acrl, PSIfull, PSIcrl, PSIdiv, Pxmat);

% Note that RBFs can't capture constant fields very well, so make sure that
% the field isn't nearly constant (i.e., the zero field) before calling.
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uu1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uu2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uu3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
hold on
quiv = quiver3(xxq,yyq,zzq,uu1,uu2,uu3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
hold off
end

