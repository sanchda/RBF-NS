function Z = getDivFree(L, X)

Xx = repmat(X(:,1),1,2*L+1);
Xy = repmat(X(:,2),1,2*L+1);
Xz = repmat(X(:,3),1,2*L+1);

[~, Dx, Dy, Dz] = dsph(L,X(:,1),X(:,2),X(:,3));

Zx = (1/sqrt(L^2+L))*(-Xy.*Dz + Xz.*Dy);
Zy = (1/sqrt(L^2+L))*(Xx.*Dz - Xz.*Dx);
Zz = (1/sqrt(L^2+L))*(-Xx.*Dy + Xy.*Dx);

Z = reshape([Zx(:) Zy(:) Zz(:)]',3*size(Zx,1), [])';

Z = reshape(Z,[size(Zx,2), 3, size(Zx,1)]);

end

