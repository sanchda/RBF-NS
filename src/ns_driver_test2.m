%% Outer loop and constants
Nvect = 10:5:80;

nu=1/10;
omega=0;

matlabpool(12);

% Empirically determined Leray projector shape parameter
% TODO: find better fit for this; especially for higher N
divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);

% Differential operator shape parameter
beta = 12;
c = (4*pi)^2;
PDE_geteps = @(N) (beta/c)*(N+1)^(9/8);

% Loop through N.  Everything has to be re-initialized between N, so no foul here.
for N = Nvect
    %% Initialize everything
    eps_Leray = divFree_geteps(N);
    eps_PDE   = PDE_geteps(N);

    % Generate RBF nodes
    X = getMDPointsNix(N);
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
    % Initialize Grady's test 2 stuff
    atm = setupT4(X);
    tf = 2*pi*atm.a/atm.u0; % seconds
    t = linspace(0,tf,23);
    disp('Test matrices initialized')

    omega=atm.omega;
    nu = atm.nu;

    % Initialize operators
    % TODO: increase speed
    [lap, projgrad, Lx, Ly, Lz, Achol, Aleray, Pxmat] = nsInitS2(X, eps_Leray, eps_PDE);
    disp('NS working matrices initialized')

    % Outer simulation loop; iterates through parameters

    
    % Timestep width (in seconds).
    h=1200;
    
    % Max time.
    maxt = 5*24*60*60; % five days
    maxc = maxt/h;
    errmat = zeros(maxc,1);
    errrel = errmat;
    
    % Get initial
    U = computeTrueSolution(atm,0);
    t=0;

    for c = 2:maxc
        [U,t] = navierstokes_test2(X, U, h, t, 1, 1, omega, lap, projgrad, Aleray, Pxmat, atm);
        Utest2 = computeTrueSolution(atm,t);

        errU = U - Utest2;
        errmag = sqrt(errU(:,1).^2 + errU(:,2).^2 + errU(:,3).^2);
        Umag = sqrt(U(:,1).^2 + U(:,2).^2 + U(:,3).^2);
        Utest2mag = sqrt(Utest2(:,1).^2 + Utest2(:,2).^2 + Utest2(:,3).^2);
        
        Umagmin=min(Umag,Utest2mag);
        mask = Umagmin > (10^0.5);
        
        errmat(c) = mean(errmag(mask));
        errrel(c) = mean(errmag(mask)./Umagmin(mask));
        
             [xx,yy,zz] = sphere(100);
     surf(xx,yy,zz,1+0*xx,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
     hold on;
     set(gcf,'color','none')
     quiver3(X(:,1),X(:,2),X(:,3), U(:,1)/     100, U(:,2)/     100, U(:,3)/     100, 0,'b-');
     quiver3(X(:,1),X(:,2),X(:,3), Utest2(:,1)/100, Utest2(:,2)/100, Utest2(:,3)/100, 0,'r-');
     hold off;
     daspect([1 1 1]);
     camorbit(140,0)
     drawnow;

        if(mod(c,50) == 0)
            sprintf('N: %d, step: %d',N,c)
        end
    end

% semilogy((1:maxc)/72,errrel)
% xlabel('time (days)','interpreter','latex','fontsize',16)
% ylabel('relative $\ell_2$ error','interpreter','latex','fontsize',16)
% title('relative $\ell_2$ error, $N=1681$','interpreter','latex','fontsize',16)


     [xx,yy,zz] = sphere(100);
     surf(xx,yy,zz,1+0*xx,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
     hold on;
     set(gcf,'color','none')
     quiver3(X(:,1),X(:,2),X(:,3), U(:,1)/     100, U(:,2)/     100, U(:,3)/     100, 0,'b-');
     quiver3(X(:,1),X(:,2),X(:,3), Utest2(:,1)/100, Utest2(:,2)/100, Utest2(:,3)/100, 0,'r-');
     hold off;
     daspect([1 1 1]);
     camorbit(140,0)
     drawnow;
    
    
    suffix     = sprintf('%d_.mat', N);
    errmat_str = sprintf('errmat_%s', suffix);
    errrel_str = sprintf('errrel_%s', suffix);  

    save(errmat_str,'errmat');
    save(errrel_str,'errrel');
    
    clear errrel;
    clear errmat;
          
end % N loop

save(eps_total,'eps_total.mat');
