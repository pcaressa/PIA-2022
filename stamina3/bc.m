%funcion BC
function BC = bc(x)
    if x <= 900
        BC = 1;
    elseif 900 < x && x <= 1800
        BC=(x/450) - (1/810000)*(x^2);
    else 
        BC=0;
    endif
    
    BC = BC * smax
