%cd 'C:\Users\david\Documents\GitHub\RBF-NS\src';

%% Initialization
cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';

%==========================================================================
%                         Parameters and Constants                        
%==========================================================================

nu = 1/10000;      % Parameter for the NS equation
omega = 1;     % Strength of coriolis force
N = 23;        % Somehow related to the number of centers.  For the
               % ME points, the number of centers is (N+1)^2.
N0 = 1;        % Highest spherical harmonic in the test
M = 1;         % how many iterations to run the simulation for
h = 1/(N+1);   % timestep
divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);
eps_Leray = divFree_geteps(N);
eps_PDE   = 1.75;
surfeps = 2;

disp('Set parameters')
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

disp('RBF nodes loaded')
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

[lap grad Lx Ly Lz Achol Afull Acrl Adiv PSIfull PSIcrl PSIdiv Pxmat] = nsInitS2(X, HGA, eps_Leray, eps_PDE);
disp('NS working matrices initialized')
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
U0 = makeGaneshTest1(N0, X, 0, nu);
%U0  = getDivFree(2,X); 
%U0  = U0(:,1:3);

disp('Initial VF generated')
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
sz = [201, 351];
szq = [35 35];
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

disp('Visualization parameters set')
%% Simulate with RBF
U = U0;

% Initial camera parameters
cpos = [-0.031 -21.392 9.115];
ctarg = [-0.031 0.122 -0.009];
cview = 4.0;

for c = 1:50
% Start with the visualization first, so the initial condition is seen.

% Note that RBFs can't capture constant fields very well, so make sure that
% the field isn't nearly constant (i.e., the zero field) before calling.
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results

clf
set(gca,'CameraPosition',cpos)
set(gca,'CameraTarget',ctarg)
set(gca,'CameraViewAngle',cview)
%caxis([-0.2,0.2])
colorbar('EastOutside')
cpos = campos;
ctarg = camtarget;
hold on
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
box off
axis off
hold off
text(-1.23,0,1.05,sprintf('t: %f',c*h),'Fontsize',12)
text(-1.23,0,0.95,sprintf('max: %f',max(max(U))),'Fontsize',12)
text(-1.23,0,0.85,sprintf('min: %f',min(min(U))),'Fontsize',12)
text(-1.23,0,-0.75,sprintf('n: %d',size(U,1)),'Fontsize',12)
text(-1.23,0,-0.85,sprintf('nu: %f',nu),'Fontsize',12)
text(-1.23,0,-0.95,sprintf('omega: %f',omega),'Fontsize',12)

F(c) = getframe(gcf);

U = navierstokes(X, U, h, 1, eps, nu, omega, N0, lap, grad, surfeps, Lx, Ly, Lz, Afull, Acrl, PSIfull, PSIcrl, PSIdiv, Pxmat);
end


movie(gcf,F,1);
%% Save movie
movie2avi(F, 'NS_576DAS_7.4.13_trial7_4.avi')

%% View corresponding Ganesh solution
U = U0;

% Initial camera parameters
cpos = [-0.031 -21.392 9.115];
ctarg = [-0.031 0.122 -0.009];
cview = 4.0;

for c = 1:50
t = h*c;

% Note that RBFs can't capture constant fields very well, so make sure that
% the field isn't nearly constant (i.e., the zero field) before calling.
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq);
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq);
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq);

% Plot the results

clf
set(gca,'CameraPosition',cpos)
set(gca,'CameraTarget',ctarg)
set(gca,'CameraViewAngle',cview)
%caxis([-0.2,0.2])
colorbar('EastOutside')
cpos = campos;
ctarg = camtarget;
hold on
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
box off
axis off
hold off
text(-1.23,0,1.05,sprintf('t: %f',c*h),'Fontsize',12)
text(-1.23,0,0.95,sprintf('max: %f',max(max(U))),'Fontsize',12)
text(-1.23,0,0.85,sprintf('min: %f',min(min(U))),'Fontsize',12)
text(-1.23,0,-0.75,sprintf('n: %d',size(U,1)),'Fontsize',12)
text(-1.23,0,-0.85,sprintf('nu: %f',nu),'Fontsize',12)
text(-1.23,0,-0.95,sprintf('omega: %f',omega),'Fontsize',12)

G(c) = getframe(gcf);


U = makeGaneshTest1(N0, X, t, nu);
end

movie(gcf,G,1)

%% Save movie
movie2avi(G(1:50), 'NS_576GAN_7.4.13_trial7.avi')
