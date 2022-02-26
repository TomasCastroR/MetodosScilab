// Practica 2

clc // limpia la consola
clear // borra el contenido de la memoria
// Primera funcion
function y = P1(x)
y = x.^7 - 7*x.^6 + 21*x.^5 - 35*x.^4 + 35*x.^3 - 21*x.^2 + 7*x - 1;
endfunction
// Segunda funcion
function y = P2(x)
y = (x - 1).^7;
endfunction
// Evaluacion de ambas funciones cerca de uno
x = linspace(1-1e-2,1+1e-2,2001);
y1 = P1(x);
y2 = P2(x);
// Grafica de las funciones
//plot(x,y1,'b');
//plot(x,y2,'r','thickness',2)
//legend(["$P1(x)$";"$P2(x)$"]);


function r = misraices(p)
    c = coeff(p,0);
    b = coeff(p,1);
    a = coeff(p,2);
    if b < 0 then
        r(1) = (-b + sqrt(b^2-4*a*c))/(2*a);
        r(2) = (2*c)/(-b + sqrt(b^2-4*a*c));
    else
        r(1) = (2*c)/(-b - sqrt(b^2-4*a*c));
        r(2) = (-b - sqrt(b^2-4*a*c))/(2*a);
    end 
endfunction

p = poly([-0.0001 10000.0 0.0001],"x","coeff");
valorEsperado = 1e-8;

raicesMias = misraices(p);
raicesScilab = roots(p); 

errorMio = abs(raicesMias(1)-valorEsperado)/valorEsperado;
errorScilab = abs(raicesScilab(2)-valorEsperado)/valorEsperado;

printf("Esperado : %e\n", valorEsperado);
printf("misraices (nuestro) : %e  (error= %e)\n", raicesMias(1), errorMio)
printf("roots (Scilab) : %e (error= %e)\n", raicesScilab(2), errorScilab);
