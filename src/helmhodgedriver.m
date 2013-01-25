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



                
    % Epsilon was fit against N
    eps = @(N) -1.37119 + 0.556389*sqrt(N);
    N   = 20;         % Somehow related to the number of centers.  For the ME points,
                      % the number of centers is (N+1)^2.
    epsilon = eps(N);    % Shape paramater for the RBF kernel


     errinfmat = zeros(70-7);
     err2mat   = zeros(70-7);
     kappamat  = zeros(70-7);
     epsmat    = zeros(70-7);

    for N = 7:30
        epsilon = eps(N);
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
        %U2 = [X(1,:).*(-X(2,:).^2 + X(3,:).^2) X(2,:).*(X(1,:)-X(3,:)).*(X(1,:) + X(3,:)), (-X(1,:).^2 + X(2,:).^2).*X(3,:)];

        % Test field 3
        %U3 = [X(:,3) -exp(X(:,1)).*X(3,:), -X(1,:) + exp(X(1,:)).*X(2,:)];

        % Run the divergence-free test
        [maxerr l2err kappa] = testHelmHodge(X, W, U1div, U1crl, HGA, epsilon);

    %    disp([maxerr l2err kappa])
        errinfmat(N-6) = maxerr;
        err2mat(N-6)   = l2err;
        kappamat(N-6)  = kappa;
        epsmat(N-6)    = epsilon;

    end