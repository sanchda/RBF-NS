% cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';
%cd 'C:\Users\david\Documents\GitHub\RBF-NS\src';
cd '/home/sanchez/RBF-NS/src'; 

eps = @(N) -0.519226 + 0.106809*(N+1);
N   = 10;           % Somehow related to the number of centers.  For the ME points,
                   % the number of centers is (N+1)^2.
 
 %Projection operator
 Q = @(x) [0 x(3) (-x(2)); (-x(3)) 0 x(1);  x(2) (-x(1)) 0];

 c=0;

 for N = 5:70
 
    epsilon = eps(N);
    disp(N);
    X = getMDPoints(N);
    W = X(:,4);
    X = X(:,1:3);

    % Rotate X through by a small angle
    t=0.5;
    theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
    for i = 1:(N+1)^2
        X(i,:) = (theta*X(i,:)')';
    end
    
    % A div-free VSH.  Check VecSphHarm.nb for details
    U1 = [X(:,1).*X(:,2) (X(:,3).^2 - X(:,2).^2) -X(:,2).*X(:,3)];
    U1 = [0*X(:,1) X(:,3) -X(:,2)];
    U1 = [X(:,3) -X(:,2) 0*X(:,1)];
    
    % Generate a divergence-free VSH from Grady's code.
    mu = 3;
    [Y Dx Dy Dz] = dsph(mu,X(:,1),X(:,2),X(:,3));
    nu = floor(mu/2);
    Dx = Dx(:,nu);
    Dy = Dy(:,nu);
    Dz = Dz(:,nu);
    U1 = 0*X;
    for i = 1:size(X,1)
        U1(i,:) = (Q(X(i,:))*[Dx(i) Dy(i) Dz(i)]')';
    end
    
    % Run the test!
    [maxerr l2err kappa] = testVecLap(X, W, U1, epsilon);
    disp([maxerr l2err])
    kappa

    c = c+1;
    errinfmat(c) = maxerr;
    err2mat(c)   = l2err;
    kappamat(c)  = kappa;
    epsmat(c)    = epsilon;
    
end
