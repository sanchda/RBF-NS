function Ucrl = projectCrlFree(U, dmat, emat, Afull, PSIcrl)
    % AUTHOR:   David Sanchez
    % DATE:     May 2013
    % MODIFIED: 7/4/2013
    
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
    Ucrl = PSIcrl*reshape(coeffs',size(coeffs,1)*size(coeffs,2),[]);
    Ucrl = reshape(Ucrl,3,[])';

end

