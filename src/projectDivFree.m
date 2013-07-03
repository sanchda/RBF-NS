function Udiv = projectDivFree(U, dmat, emat, Afull, PSIdiv)
    % AUTHOR:   David Sanchez
    % DATE:     May 2013
    % MODIFIED: 7/3/2013
    
    % Perform RBF-based projection onto the divergence-free component of a
    % vector field on a boundary-free 2-manifold.
    %
    %-------------------------------ARGUMENTS------------------------------
    %
    % U is the vector-field, a matrix with 3 columns and N rows, where N is
    % the number of RBF centers on the manifold.  Implicitly, there should
    % be a matrix X such that f(X(i,:)) = U(i,:) and X are the coordinates
    % of the centers in R3.
    % ASSERT:  U is assumed to be in the tangent space of the manifold!
    % Non-tangential vector-fields are unsupported and unchecked for.
    %
    % dmat, emat are matrices holding the local coordinate frames for the
    % values in X.
    % i.e., [dmat emat] are the local coordinates of X.
    % 
    % Afull is the matrix of vector RBF coefficients (A^(2) in
    % Fuselier&Wright's notation) containing the coefficients for both the
    % divergence- and curl-free components.
    %
    % PSIdiv is the matrix of the divergence-free matrix RBF kernel.
    % i.e., each i-j 3x3 block of PSIdiv is the result of applying the
    % divergence-free component of the function \Psi_{div} in the
    % literature to X(i,:) and X(j,:).  This matrix is built in
    % navierstokes.m
    %
    % Udiv is the divergence-free component of U, given in the same
    % coordinates as U.
    %
    %-------------------------------BEHAVIOR-------------------------------
    %
    % For mathematical details, refer to the Fuselier&Wright (and Narcowich
    % and Ward) papers on surface-divergence free RBF interpolants.
    % http://math.boisestate.edu/~wright/research/
    %
    %
    
   %% Set up Ax = b (determine b) 
    % Get coefficients
    for i = 1:size(U,1)
        gamdel((2*i-1):(2*i)) = [dmat(i,:);emat(i,:)]*(U(i,:).');
    end

    % Solve for x
    albet = Afull\gamdel';

    coeffs = repmat(albet(1:2:size(albet,1)),1,3).*dmat + ...
             repmat(albet(2:2:size(albet,1)),1,3).*emat;

    % Finally, recover the div-free part.
    Udiv = PSIdiv*reshape(coeffs',size(coeffs,1)*size(coeffs,2),[]);
    Udiv = reshape(Udiv,3,[])';

end

