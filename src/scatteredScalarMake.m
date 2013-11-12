function Wnew = scatteredScalarMake(X, Y, Z, W, Xnew, Ynew, Znew, eps)

% Tentatively, truncates values outside of the unit sphere!!

A = (repmat(X,[1 size(X,1)]) - repmat(X',[size(X,1) 1])).^2;
A = A + (repmat(Y,[1 size(X,1)]) - repmat(Y',[size(X,1) 1])).^2;
A = A + (repmat(Z,[1 size(X,1)]) - repmat(Z',[size(X,1) 1])).^2;
A = exp(-eps*eps*A);
coeffs = (A')\W;

Wnew = 0*Znew;

for i=1:size(Xnew,1)
    for j = 1:size(Ynew,1)
        p = [Xnew(i,j) Ynew(i,j) Znew(i,j)];
        Anew = (X - p(1)).^2 + (Y - p(2)).^2 + (Z - p(3)).^2;
        Anew = exp(-eps*eps*Anew);
        Wnew(i,j) = dot(coeffs,Anew);
    end
end


end

