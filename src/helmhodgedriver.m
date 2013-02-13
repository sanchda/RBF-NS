    cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';
    % GA RBF
    HGA =  @(x,y,eps) [exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
                    -1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(1)+(-1).*y(1)).^2),4.*exp( ...
                    1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).* ...
                    y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(2)+(-1).*y(2)),4.*exp(1).^((-1).* ...
                    eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).* ...
                    eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3));4.*exp(1).^((-1).*eps.^2.*(( ...
                    x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(1)+( ...
                    -1).*y(1)).*(x(2)+(-1).*y(2)),exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+( ...
                    x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(2)+( ...
                    -1).*y(2)).^2),4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1) ...
                    .*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)); ...
                    4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
                    -1).*y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3)),4.*exp(1).^(( ...
                    -1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)) ...
                    .*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)),exp(1).^((-1).*eps.^2.*(( ...
                    x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).* ...
                    eps.^2+4.*eps.^4.*(x(3)+(-1).*y(3)).^2)];

    j=0;
    
    % The function Q(x) generates the matrix which will project a vector in R3
    % onto the tangent space of the sphere at x \in S2
    Q = @(x) [0 x(3) (-x(2)); (-x(3)) 0 x(1);  x(2) (-x(1)) 0];
    % P(x) projects into the normal space
    P = @(x) eye(3) - [x(1);x(2);x(3)]*[x(1) x(2) x(3)];
    
    for N = 5:5:35
        geteps = @(n) -0.519226 + 0.106809*(n+1);
        epsilon = geteps(N);
        X = getMEPoints(N);
        W = X(:,4);
        X = X(:,1:3);

        % Rotate X through by a small angle
        t=0.5;
        theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
        for i = 1:(N+1)^2
            X(i,:) = (theta*X(i,:)')';
        end

        % Div-free VFs can be given by curls of scalars.  Just use Qx*grad(u) =
        % curl_T(u) on the sphere.

        % Curl-free VFs

        % Test field 1, div-free component
        U1div = [(-X(:,2) + X(:,3)) (X(:,1) - X(:,3)) (-X(:,1) + X(:,2))];
        % Test field 1, curl-free component
        U1crl = [1-X(:,1).*(X(:,1) + X(:,2) + X(:,3)),1-X(:,2).*(X(:,1) + X(:,2) + X(:,3)),1-X(:,3).*(X(:,1) + X(:,2) + X(:,3))];

        % Test field 2
        U2crl = [(-2).*exp(1).^X(:,3).*((-1)+X(:,1).^2).*cos(2.*X(:,1))+((-2)+2.*X(:,1).^2+(-1).*exp(1).^X(:,3).*X(:,1).*X(:,3)).*sin(2.*X(:,1)),(-1).*X(:,2).*(2.*exp(1).^X(:,3).*X(:,1).*cos(2.*X(:,1))+((-2).*X(:,1)+exp(1).^X(:,3).*X(:,3)).*sin(2.*X(:,1))),(-2).*exp(1).^X(:,3).*X(:,1).*X(:,3).*cos(2.*X(:,1))+(2.*X(:,1).*X(:,3)+(-1).*exp(1).^X(:,3).*((-1)+X(:,3).^2)).*sin(2.*X(:,1))];
        U2div = [(-1).*exp(1).^X(:,3).*X(:,2).*sin(2.*X(:,1)),(-2).*exp(1).^X(:,3).*X(:,3).*cos(2.*X(:,1))+(exp(1).^X(:,3).*X(:,1)+2.*X(:,3)).*sin(2.*X(:,1)),2.*X(:,2).*(exp(1).^X(:,3).*cos(2.*X(:,1))+(-1).*sin(2.*X(:,1)))];
        
        % Test field 3
        [Y,Dx,Dy,Dz] = dsph(10,X(:,1),X(:,2),X(:,3));
        U3div = 0*X;
        U3crl = 0*X;
        Dx = Dx(:,3);
        Dy = Dy(:,3);
        Dz = Dz(:,3);
        for i = 1:size(X,1)
           U3div(i,:) = (Q(X(i,:))*[Dx(i,:);Dy(i,:);Dz(i,:)])';
           U3crl(i,:) = (P(X(i,:))*[Dx(i,:);Dy(i,:);Dz(i,:)])';
        end
        

        % Run the divergence-free test
        [maxerr l2err kappa] = testHelmHodge(X, W, U3div, U3crl, HGA, epsilon);
        disp([maxerr l2err])
        j = j+1;
        errinfmat(j) = maxerr;
        err2mat(j)   = l2err;
        kappamat(j)  = kappa;
        epsmat(j)    = epsilon;

    end