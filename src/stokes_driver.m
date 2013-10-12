%==========================================================================
%                         Parameters and Constants                        
%==========================================================================
for N = 10:10:50
	disp('Initializing new N')
	disp(N)

	nu = 1/10;      % Parameter for the NS equation
	omega = 0;     % Strength of coriolis force

	% Generated parameters
	divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);
	eps_Leray = divFree_geteps(N);
	beta = 12;
	c = (4*pi)^2;
	eps_PDE = (beta/c)*(N+1)^(9/8);
	eps_Leray = eps_Leray;
	surfeps = 2;
	
	disp('Parameters Set')
	%==========================================================================
	%                            Generate RBF nodes                        
	%==========================================================================

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


	%==========================================================================
	%                            Hessian for matrix RBF                        
	%==========================================================================
	% I am SO sorry for this.  For ease of use, it was generated using a
	% Mathematica script that converts Mathematica output form to Matlab input.
	%
	% TODO:  generate using Matlab's symbolic toolkit
	
	% Moved to own file, rbf_HGA

	% This is an amazingly wasteful way of defining the operators used
	% throughout the NS code, but it has the advantage of providing a
	% manifold-agnostic interface.
	% TODO: 
	% * streamline and clean up
	% * determine if any of these have a sparse structure, make sparse
	% * vectorize initialization code

	[lap projgrad Lx Ly Lz Achol Aleray Pxmat] = nsInitS2(X, eps_Leray, eps_PDE);
	disp('NS working matrices initialized')
	%==========================================================================
	%                            Generate initial VF                       
	%==========================================================================
	h = 1/200;
	for N0 = [1 2 3 4 5 10]
	  disp('New N0')
	  disp(N0)

	  U0 = makeGaneshTest1(N0, X, h, nu);
	


	  %% Simulate with RBF
	  h=1/200;
	  U0 = makeGaneshTest1(N0, X, h, nu);
	  U = U0;
	  Uganesh = U0;
	  t=h;

	  maxc=2000;
	  errmat = zeros(maxc,1);
	  errmax = errmat;
	  tgan   = zeros(maxc-1,1);
	  tnav   = tgan;

	  for c = 2:maxc

	    tic;
	    Uganesh = makeGaneshTest1(N0, X, t+h, nu);
	    tgan(c-1) = toc;
	    tic;
	    [U,t] = navierstokes(X, U, h, t, 1, nu, omega, N0, lap, projgrad, Aleray, Pxmat);
	    tnav(c-1) = toc;
	    errmean = U - Uganesh;
	    errmean = sqrt(errmean(:,1).^2 + errmean(:,2).^2 + errmean(:,3).^2);
	    errmax(c) = max(errmean);
	    errmat(c) = mean(errmean);
	    if(mod(c,100) == 0)
		 disp(c)
	    end
	  end

	  suffix     = sprintf('%d_%d.mat', N, N0);
	  tgan_str   = sprintf('tgan_%s', suffix);
	  tnav_str   = sprintf('tnav_%s', suffix);
	  errmat_str = sprintf('errmat_%s', suffix);
	  errmax_str = sprintf('errmax_%s', suffix);  

	  save(tgan_str,'tgan');
	  save(tnav_str,'tnav');
	  save(errmat_str,'errmat');
	  save(errmax_str,'errmax');

	  disp 'Error for current round:'
	  errmat(end)
	  clear tgan;
	  clear tnav;
	  clear errmax;
	  clear errmat;

	end
end
