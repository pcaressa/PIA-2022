function G = grad3D(k, m, dx, n, dy, o, dz)
% Returns a three-dimensional mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%                o : Number of cells along z-axis
%               dz : Step size along z-axis

    Im = sparse(m + 2, m);
    Im(2:(m + 2) - 1, :) = speye(m, m);
    
    Gx = grad(k, m, dx);
    
    In = sparse(n + 2, n);
    In(2:(n + 2) - 1, :) = speye(n, n);
    
    Gy = grad(k, n, dy);
    
    Io = sparse(o + 2, o);
    Io(2:(o + 2) - 1, :) = speye(o, o);
    
    Gz = grad(k, o, dz);
    
    Sx = kron(kron(Io', In'), Gx);
    Sy = kron(kron(Io', Gy), Im');
    Sz = kron(kron(Gz, In'), Im');
    
    G = [Sx; Sy; Sz];
end
