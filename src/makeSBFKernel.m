function Asbf = makeSBFKernel(X, PSI, d, e)
% Generates the matrix A in order to solve for the coefficients of f_i on
% the tangent space of the sphere.

% X         Nx3 array of cell centers, which are vectors on S2 in R3
% PSI       Handle to a function computing the divergence-free interpolant
% d,e       Handles to zonal and meridional bases


N = size(X,1);
Asbf = zeros(2*N,2*N);

% matrix with all of the D on top and all of the E on bottom
% configure for simple indexing

dmat = d(X(:,:));
emat = e(X(:,:));

leftmat = reshape([dmat(:) emat(:)]',2*size(dmat,1),[]);
leftmat = reshape(leftmat',3,2,[]);
leftmat = permute(leftmat, [2 1 3]);


% Can either use nested for loops:
%
% for i = 1:size(X,1)
%   for j = 1:size(X,1)
%       Asbf( (2*i-1):(2*i),(2*j-1):(2*i) ) = [dmat(i) emat(i)]*PSI(X(j,:),X(i,:))*[dmat(j)';emat(j)']'
%    end
% end
%
% Or, use arrayfun to collapse one of the nested loops
%
for i =1:size(X,1)
    Acol = arrayfun(@(j) ...
    	leftmat(:,:,j)*PSI(X(j,:), X(i,:))*leftmat(:,:,i)', (1:size(X,1))','UniformOutput', 0);
    Asbf(:,(2*i-1):(2*i)) = cell2mat(Acol);

end

end

