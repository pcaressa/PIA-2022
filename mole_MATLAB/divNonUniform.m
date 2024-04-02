function D = divNonUniform(k, ticks)
% Returns a m+2 by m+1 one-dimensional non-uniform mimetic divergence
% operator
%
% Parameters:
%                k : Order of accuracy
%                ticks : Edges' ticks e.g. [0 0.1 0.15 0.2 0.3 0.4 0.45]

    % Get uniform operator without scaling
    D = div(k, length(ticks)-1, 1);
    
    [m, ~] = size(D);
    
    % Compute the Jacobian using the uniform operator and the ticks
    if size(ticks, 1) == 1
        J = spdiags((D*ticks').^-1, 0, m, m);
    else
        J = spdiags((D*ticks).^-1, 0, m, m);
    end
    
    % This is the non-uniform operator
    D = J*D;
end
