% cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';
cd 'C:\Users\david\Documents\GitHub\RBF-NS\src';
            
eps = @(N) -1.37119 + 0.556389*sqrt(N);
epsilon = 2.05;    % Shape paramater for the RBF kernel
N   = 10;           % Somehow related to the number of centers.  For the ME points,
                   % the number of centers is (N+1)^2.

                 
 errinfmat = zeros(70-7);
 err2mat   = zeros(70-7);
 kappamat  = zeros(70-7);
 epsmat    = zeros(70-7);
 
 %for N = 7:70
 
    epsilon = eps(N);
    disp(N);
    X = getMEPoints(N);
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
    
    % Run the test!
    [maxerr l2err kappa] = testVecLap(X, W, U1, epsilon);

    disp([maxerr l2err kappa])
    errinfmat(N-6) = maxerr;
    err2mat(N-6)   = l2err;
    kappamat(N-6)  = kappa;
    epsmat(N-6)    = epsilon;
    
%end