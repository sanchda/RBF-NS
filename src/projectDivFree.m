function Udiv = projectDivFree(U, Aleray)
    % AUTHOR:   David Sanchez
    % DATE:     May 2013
    % MODIFIED: 7/3/2013
    
    % Perform RBF-based projection onto the divergence-free component of a
    % vector field on a boundary-free 2-manifold.
    %
    % OBSERVATION:
    % This function should be idempotent, up to numerical accuracy (which
    % will be around machine precision for n>=12^2 on the sphere).
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
    % Aleray is a specially-crafted projection operator, taken by observing
    % that if Uc, Us are vectors of the same "points", but in Cartesian and
    % spherical coordinates resp and Acs, Asc take cartesian to sph and sph
    % to cart, resp., Asc L Acs U_c = U_div_c.
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
     Udiv = reshape((Aleray*reshape(U',[],1)),3,[])';
end

