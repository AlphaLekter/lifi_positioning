function [f, J] = sphere_intersect(x, c1, r1, c2, r2, c3, r3, c4, r4, c5, r5, number_of_sphere)
% Ora possiamo definire una funzione che calcola la differenza tra la distanza 
% tra due punti (i.e., il centro delle rispettive sfere) e i loro raggi, rappresentare 
% ciò come un vettore e calcolare una norma del quadrato di tale vettore. 
% Questa funzione sarà poi utilizzata come funzione di costo per l'algoritmo
% di Levenberg-Marquardt.


if number_of_sphere == 5
    %calcolo le distanze
    d1 = norm(x - c1) - r1;
    d2 = norm(x - c2) - r2;
    d3 = norm(x - c3) - r3;
    d4 = norm(x - c4) - r4;
    d5 = norm(x - c5) - r5;
    
    %calcolo la funzione di costo
    f = [d1; d2; d3; d4; d5];
    
    %calcolo la matrice Jacobiana
    J = [(x-c1)'/norm(x-c1), -(x-c2)'/norm(x-c2), -(x-c3)'/norm(x-c3), ...
        (x-c4)'/norm(x-c4), (x-c5)'/norm(x-c5)];
elseif number_of_sphere == 4
        %calcolo le distanze
    d1 = norm(x - c1) - r1;
    d2 = norm(x - c2) - r2;
    d3 = norm(x - c3) - r3;
    d4 = norm(x - c4) - r4;
   
    
    %calcolo la funzione di costo
    f = [d1; d2; d3; d4];
    
    %calcolo la matrice Jacobiana
    J = [(x-c1)'/norm(x-c1), -(x-c2)'/norm(x-c2), -(x-c3)'/norm(x-c3), ...
        (x-c4)'/norm(x-c4)];
elseif number_of_sphere == 3
        %calcolo le distanze
    d1 = norm(x - c1) - r1;
    d2 = norm(x - c2) - r2;
    d3 = norm(x - c3) - r3;
    
    %calcolo la funzione di costo
    f = [d1; d2; d3];
    
    %calcolo la matrice Jacobiana
    J = [(x-c1)'/norm(x-c1), -(x-c2)'/norm(x-c2), -(x-c3)'/norm(x-c3)];    
else
    f = NaN;
    J = NaN;
end

end