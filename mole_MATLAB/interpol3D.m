function I = interpol3D(m, n, o, c1, c2, c3)
% Returns a three-dimensional interpolator of 2nd-order
%                m : Number of cells along x-axis
%                n : Number of cells along y-axis
%                o : Number of cells along z-axis
%               c1 : Left interpolation coeff.
%               c2 : Bottom interpolation coeff.
%               c3 : Front interpolation coeff.

    Im = sparse(m + 2, m);
    Im(2:(m + 2) - 1, :) = speye(m, m);
    
    Ix = interpol(m, c1);
    
    In = sparse(n + 2, n);
    In(2:(n + 2) - 1, :) = speye(n, n);
    
    Iy = interpol(n, c2);
    
    Io = sparse(o + 2, o);
    Io(2:(o + 2) - 1, :) = speye(o, o);
    
    Iz = interpol(o, c3);
    
    Sx = kron(kron(Io', In'), Ix);
    Sy = kron(kron(Io', Iy), Im');
    Sz = kron(kron(Iz, In'), Im');
    
    I = [Sx; Sy; Sz];
end
