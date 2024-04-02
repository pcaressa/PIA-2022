function S = sidedNodal(m, dx, type)
% Returns a m+1 by m+1 one-dimensional sided approximation for uniformly
% spaced data points. This function is handy for advective terms.
%
% Parameters:
%                m : Number of cells
%               dx : Step size
%             type : 'backward', 'forward' or 'centered'

    switch type
        case 'backward'
            S = spdiags([-ones(m+1, 1) ones(m+1, 1)], [-1 0], m+1, m+1);
            S(1, end-1) = -1;
            S = S/dx;
        case 'forward'
            S = spdiags([-ones(m+1, 1) ones(m+1, 1)], [0 1], m+1, m+1);
            S(end, 2) = 1;
            S = S/dx;
        otherwise % 'centered'
            S = spdiags([-ones(m+1, 1) zeros(m+1, 1) ones(m+1, 1)], [-1 0 1], m+1, m+1);
            S(1, end-1) = -1;
            S(end, 2) = 1;
            S = S/(2*dx);
    end
end
