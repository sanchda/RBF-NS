function Z = getDivFree(L, X)
% AUTHOR:   David Sanchez
% DATE:     August 2012
% MODIFIED: 7/5/2013

% Produces the divergence-free vector spherical harmonics of order L,
% evaluated at the points in X, assumed to be a 3-column matrix of floats
% corresponding to 3-tuples, (x,y,z) with x^2+y^2+z^2==1 (on the sphere).

% Output is formatted for compatibility with Grady Wright's dsph code,
% which produces the scalar spherical harmonics and their derivatives.
%
% The first three columns of Z are the order L div-free VSH of index -L
% The next three columns forms the index -(L-1) VSH
% The index 0 VSH is at Z(:,(3*L+1):(3*L+4)).
%
% Some mathematical properties of the output:
% Each div-free VSH is divergence-free on the sphere (duh) and an
% eigenfunction of the vector Laplacian on the sphere.

Xx = repmat(X(:,1),1,2*L+1);
Xy = repmat(X(:,2),1,2*L+1);
Xz = repmat(X(:,3),1,2*L+1);

[~, Dx, Dy, Dz] = dsph(L,X(:,1),X(:,2),X(:,3));

Zx = ((1/sqrt(L^2+L))*(-Xy.*Dz + Xz.*Dy))';
Zy = ((1/sqrt(L^2+L))*(Xx.*Dz - Xz.*Dx))';
Zz = ((1/sqrt(L^2+L))*(-Xx.*Dy + Xy.*Dx))';

Z = reshape([Zx(:) Zy(:) Zz(:)]',3*size(Zx,1),[])';

end

