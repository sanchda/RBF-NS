function [Asbf, Asbf_test] = makeSBFKernel(X, d ,e, eps_Leray, type)
% Generates the matrix A in order to solve for the coefficients of f_i on
% the tangent space of the sphere.

% X         Nx3 array of cell centers, which are vectors on S2 in R3
% PSI       Handle to a function computing the divergence-free interpolant
% d,e       Handles to zonal and meridional bases


N = size(X,1);

% We want Asbf to be a 2N x 2N array, but here we use parfor to
% parallelize, requiring sliced arrays for good performance.  Our desire is
% for each 2-unit horizontal strip of Asbf to be indexed separately.
% UPDATE: these matrices are not sliced, but merely indexed, due to the
% arrayfun.

% If not using parfor:
% Asbf = zeros(2*N, 2*N)


Asbf = zeros(N,2,2*N);
Asbf = num2cell(Asbf, [2 3]);

% Still not quite right, so reshape.
% TODO: build this in a smarter way to avoid parfor
parfor i=1:N
    Asbf{i}=reshape(Asbf{i},2,[]);
end

% matrix with all of the D on top and all of the E on bottom
% configure for simple indexing

dmat = d(X(:,:));
emat = e(X(:,:));

leftmat = reshape([dmat(:) emat(:)]',2*size(dmat,1),[]);
leftmat = reshape(leftmat',3,2,[]); % stop here if not using parfor

leftmat = permute(leftmat, [2 1 3]);

% Convert to cell array, with each cell containing dimensions 1 and 2 of
% previous array.  This defines a "sliced array", which is much faster in
% parfor.
leftmat = num2cell(leftmat, [1 2]);


switch type
    case 'full'
        parfor i =1:size(X,1)
            Asbf{i} = cell2mat(arrayfun(@(j) leftmat{j}*PSI(X(j,:), X(i,:), eps_Leray)*leftmat{i}', (1:size(X,1))','UniformOutput', 0));
        end
    case 'curl'
        parfor i =1:size(X,1)
            Asbf{i} = cell2mat(arrayfun(@(j) leftmat{j}*PSIcurl(X(j,:), X(i,:), eps_Leray)*leftmat{i}', (1:size(X,1))','UniformOutput', 0));
        end
    case 'div'
        parfor i =1:size(X,1)
            Asbf{i} = cell2mat(arrayfun(@(j) leftmat{j}*PSIdiv(X(j,:), X(i,:), eps_Leray)*leftmat{i}', (1:size(X,1))','UniformOutput', 0));
        end
    otherwise
        disp 'Something went horribly wrong'
end



% switch type
%     case 'full'
%         parfor i =1:size(X,1)
%             Acol = arrayfun(@(j) squeeze(leftmat(j,:,:))*PSI(X(j,:), X(i,:), eps_Leray)*squeeze(leftmat(i,:,:))', (1:size(X,1))','UniformOutput', 0);
%             Asbf{i} = cell2mat(Acol);
%         end
%     case 'curl'
%         parfor i =1:size(X,1)
%             Acol = arrayfun(@(j) squeeze(leftmat(j,:,:))*PSIdiv(X(j,:), X(i,:), eps_Leray)*squeeze(leftmat(i,:,:))', (1:size(X,1))','UniformOutput', 0);
%             Asbf{i} = cell2mat(Acol);
%         end
%     case 'div'
%         parfor i =1:size(X,1)
%             Acol = arrayfun(@(j) squeeze(leftmat(j,:,:))*PSIcrl(X(j,:), X(i,:), eps_Leray)*squeeze(leftmat(i,:,:))', (1:size(X,1))','UniformOutput', 0);
%             Asbf{i} = cell2mat(Acol);
%         end
%     otherwise
%         disp 'Something went horribly wrong'
% end
% 



% Asbf comes out interleaved and misshapen, so it should be fixed.
% TODO: try doing this in one reshape
Asbf=cell2mat(Asbf);

Asbf_buff=reshape(Asbf,2*N,[]);
Asbf=zeros(2*N,2*N);

Asbf(:,1:2:end)=Asbf_buff(:,1:N);
Asbf(:,2:2:end)=Asbf_buff(:,(N+1):end);

end


% The function Q(x) generates the matrix which will project a vector
% in R3 onto the tangent space of the sphere at x in S2
function y = Q(x)
y = [0*x(:,1) x(:,3) (-x(:,2));...
     (-x(:,3)) 0*x(:,1) x(:,1);...
     x(:,2) (-x(:,1)) 0*x(:,3)];
end

% explicit indexing in case numel(x)>3 for some reason
function y = P(x)
  y = eye(3) - [x(:,1);x(:,2);x(:,3)]*[x(:,1) x(:,2) x(:,3)];
end
    
% Define the matrix-valued kernel
function y = PSIdiv(x, y, eps_Leray)
  y = (Q(x)')*(-rbf_HGA(x,y,eps_Leray))*Q(y);
end

function y = PSIcrl(x, y, eps_Leray)
  y = (P(x)')*(-rbf_HGA(x,y,eps_Leray))*P(y);
end

function y = PSI(x, y, eps_Leray)
  y = PSIdiv(x, y, eps_Leray) + PSIcrl(x, y, eps_Leray);
end


function Asbf = makeSBFKernel_test(X, d, e, eps_Leray)
N = size(X,1);

Asbf = zeros(N,2,2*N);
Asbf = num2cell(Asbf, [2 3]);

% Still not quite right, so reshape.
% TODO: build this in a smarter way to avoid parfor
parfor i=1:N
    Asbf{i}=reshape(Asbf{i},2,[]);
end

% matrix with all of the D on top and all of the E on bottom
% configure for simple indexing

dmat = d(X(:,:));
emat = e(X(:,:));

leftmat2 = reshape([dmat(:) emat(:)]',2*size(dmat,1),[]);
leftmat2 = reshape(leftmat2',3,2,[]); % stop here if not using parfor

leftmat2 = permute(leftmat2, [2 1 3]);

% Convert to cell array, with each cell containing dimensions 1 and 2 of
% previous array.  This defines a "sliced array", which is much faster in
% parfor.
leftmat2 = num2cell(leftmat2, [1 2]);

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

% parfor i =1:size(X,1)
%     Acol = arrayfun(@(j) leftmat(:,:,j)*PSI(X(j,:), X(i,:))*leftmat(:,:,i)', (1:size(X,1))','UniformOutput', 0);
%     Asbf(:,(2*i-1):(2*i)) = cell2mat(Acol);
% end
% 
% end


parfor i =1:size(X,1)
    Asbf{i} = cell2mat(arrayfun(@(j) leftmat2{j}*PSI(X(j,:), X(i,:), eps_Leray)*leftmat2{i}', (1:size(X,1))','UniformOutput', 0));
end

% Asbf comes out interleaved and misshapen, so it should be fixed.
% TODO: try doing this in one reshape
Asbf=cell2mat(Asbf);

Asbf_buff=reshape(Asbf,2*N,[]);
Asbf=zeros(2*N,2*N);

Asbf(:,1:2:end)=Asbf_buff(:,1:N);
Asbf(:,2:2:end)=Asbf_buff(:,(N+1):end);

end
