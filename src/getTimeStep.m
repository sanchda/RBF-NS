function h = getTimeStep(U, minarc)
% An adaptive time-stepping strategy, based loosely on the CFL condition.

maxu = max(sqrt(U(:,1).^2 + U(:,2).^2 + U(:,3).^2));

h = min(1/sqrt(size(U,1)),0.5*minarc/maxu);

end

