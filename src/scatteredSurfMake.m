function [Xnew, Ynew, Znew] = scatteredSurfMake(X, Y, Z, eps)
% X, Y, Z should be vectors of the same size, with a functional
% relationship between Z and X x Y.

% Tentatively, truncates values outside of the unit sphere!!

A = (repmat(X,[1 size(X,1)]) - repmat(X',[size(X,1) 1])).^2;
A = A + (repmat(Y,[1 size(X,1)]) - repmat(Y',[size(X,1) 1])).^2;
A = exp(-eps*eps*A);
coeffs = (A')\Z;

[Xnew Ynew] = meshgrid(-1:2/(size(X,1)-1):1);

Znew = 0*Xnew;

for i=1:size(Xnew,1)
    for j = 1:size(Ynew,1)
        p = [Xnew(i,i) Ynew(j,j)];
        Anew = (X - p(1)).^2 + (Y - p(2)).^2;
        Anew = exp(-eps*eps*Anew);
        Znew(i,j) = dot(coeffs,Anew);
        if Xnew(i,i)^2 + Ynew(j,j)^2 > 1
            Znew(i,j) = 0;
        end
    end
end


end

