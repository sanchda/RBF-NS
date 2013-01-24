function Adiv = makeAdiv(X, PSI, d, e)
% Generates the matrix A in order to solve for the coefficients of f_i on
% the tangent space of the sphere.

% X         Nx3 array of cell centers, which are vectors on S2 in R3
% PSI       Handle to a function computing the divergence-free interpolant
% d,e       Handles to zonal and meridional bases


N = size(X,1);
Adiv = zeros(2*N,2*N);

for i=1:N
    for j=1:N
        Adiv((2*i-1):(2*i),(2*j-1):(2*j))=[e(X(i,:)); d(X(i,:))]*PSI(X(i,:),X(j,:))*[e(X(j,:)); d(X(j,:))]';
    end
end

end

