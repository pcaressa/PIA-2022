% Condizione iniziale: disco
x = i*dt;
for si = sx0:sx1
    for sj = sy0:sy1
        if (si-centro_i)^2 + (sj-centro_j)^2 <= raggio2
            % se cond1 then val1 elif cond2 then val2 else val3
            if x <= 900
                bc = smax;
            elseif 900 < x && x <= 1800
                bc = smax *((x/450) - (1/810000)*(x^2));
            else 
                bc = 0;
            endif
            S(si,sj,sz1) = bc;
            % S(si, sj, sz1) = bc(i*dt);   % b sul bordo
        endif
    endfor
endfor

%S(sx0:sx1, sy0:sy1, sz0:sz1) = b;   % b sul bordo
S(i0:i1, j0:j1, p0:p1) = 0;         % 0 nel necrotico
