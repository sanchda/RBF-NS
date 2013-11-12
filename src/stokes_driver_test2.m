%% Outer loop and constants
Nvect = 10:10:80;
N0vect = [1 2 3 5 10 20]; 
hvect  = [1/50 1/100 1/200 1/400 1/800];

nu=1/10;
omega=0;

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
    X = getMEPointsNix(N);
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
    
    % Initialize operators
    % TODO: increase speed
    [lap projgrad Lx Ly Lz Achol Aleray Pxmat] = nsInitS2(X, eps_Leray, eps_PDE);
    disp('NS working matrices initialized')
    

    % Initialize Grady's stuff
    atm = setupT4(X);
    tf = 2*pi*atm.a/atm.u0;
    t = linspace(0,tf,23);


    % Outer simulation loop; iterates through parameters
    for N0 = N0vect
        for h = hvect
          U0 = 
          U = U0;
          Uganesh = U0;
          t=h;

          % Max timesteps.  Make sure the final nondimensional time is 2
          maxc=2*(1/h);
          errmat = zeros(maxc,1);
          errmax = errmat;

          % Simulation loop
          for c = 2:maxc

            Uganesh = makeGaneshTest1(N0, X, t+h, nu);

            [U,t] = navierstokes(X, U, h, t, 1, nu, omega, N0, lap, projgrad, Aleray, Pxmat);

            errmean = U - Uganesh;
            errmean = sqrt(errmean(:,1).^2 + errmean(:,2).^2 + errmean(:,3).^2);
            errmax(c) = max(errmean);
            errmat(c) = mean(errmean);
                if(mod(c,100) == 0)
                    disp(c)
                end
          end

          suffix     = sprintf('%d_%d_%d.mat', N, N0, h);
          errmat_str = sprintf('errmat_%s', suffix);
          errmax_str = sprintf('errmax_%s', suffix);  

          save(errmat_str,'errmat');
          save(errmax_str,'errmax');


          clear tgan;
          clear tnav;
          clear errmax;
          clear errmat;
          
        end % h loop
    end % N0 loop
end % N loop

save(eps_total,'eps_total.mat');
