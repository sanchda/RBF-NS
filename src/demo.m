%
% This test was first proposed by Williamson et. al. (1992) for the shallow
% water wave equations and is called Test Case 4.  You can also find a
% description (in Cartesian coordinates) in Flyer and Wright (2009), see
% Unsteady Flow Test Case section.  I have modified it to work just for the
% incompressible Navier-Stokes equations.
%
% Note that all parameters have been set to make sense for the sphere.
% **This means that every time you compute a gradient it should be divided
% by the radius of the earth (which is stored in atm.a after setting up the
% code).**
%
% All time is measured in seconds.  I would try timesteps on the order of
% 20 minutes (1200 seconds).  
%
% Running the test to time t=5*24*60*60 (5 days) is what the test calls for.
% For it to translate around the sphere it takes about (23.1665 days).

% Load in the nodes you want to use
x = load('md042.01849'); x = x(:,1:3);

% Setup some parameters 
atm = setupT4(x);

% Code showing the true solution for one whole revolution.
tf = 2*pi*atm.a/atm.u0; % seconds
t = linspace(0,tf,23);
for j=1:length(t)
    uc = computeTrueSolution(atm,t(j));
    % Compute the force (this should be added onto the RHS of the NS equations)
    F = computeForcing(atm,t(j));

    % Plot the results
    [xx,yy,zz] = sphere(100);
    surf(xx,yy,zz,1+0*xx,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    hold on;
    quiver3(x(:,1),x(:,2),x(:,3),uc(:,1),uc(:,2),uc(:,3),'b-');
    quiver3(x(:,1),x(:,2),x(:,3),F(:,1),F(:,2),F(:,3),'r-');
    hold off;
    daspect([1 1 1]);
    drawnow;
    pause(0.2);
end
