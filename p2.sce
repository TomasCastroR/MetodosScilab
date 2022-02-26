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

// Ejercicio 1
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


// Ejemplo error relativo diferencias finitas
clc // limpia la consola
clear // borra el contenido de la memoria
xdel(winsid()) // cierra ventanas graficas
// Definicion de la funcion
function y = f(x)
y = x.*x;
endfunction
// Calculo de la derivada utilizando diferencias finitas
function y = dfa(f,x,h)
y = (f(x+h) - f(x))./h;
endfunction
x = 1; // Punto donde vamos a evaluar la derivada
ih = (0:16);
h = (10.^-ih); // Vector con los valores de h
df_approx = dfa(f,x,h); // Evaluacion de la derivada por diferencias finitas
df_scilab = numderivative(f,x,[],order=1); // Derivada obtenida por numderivative
df_true = 2; // Valor verdadero de la derivada en x = 1

// Errores absolutos y relativos
err_abs = abs(df_approx - df_true);
err_rel = err_abs./abs(df_true);
err_abs_sci = abs(df_scilab - df_true);
err_rel_sci = err_abs_sci/abs(df_true);

// Grafica
plot(ih,log10(err_rel),'b*-'); // Gr´afica en escala logar´ıtmica en el eje y
title('Error relativo utilizando diferencias finitas');
xlabel('i');
ylabel('$log {10} (Err Rel)$');
plot(ih,log10(err_rel_sci*ones(length(ih),1)),'r-');

// Impresion de resultados en pantalla
tablevalue = [ih',h',df_true*ones(length(h),1),df_approx',err_abs',err_rel'];
mprintf('%s\n',strcat(repmat('-',1,120)));
mprintf('%4s %10s \t%14s %20s \t     %16s %16s\n',...
    'i', 'h','Der. exact','Der approx','Abs. error','Rel. error');
mprintf('%s\n',strcat(repmat('-',1,120)));
mprintf('%4d %8.1e \t%9.6e %18.10e %14.5e %14.5e\n',tablevalue);
mprintf('%s\n',strcat(repmat('-',1,120)));
mprintf('%4.1s %8s \t%9.6e %18.10e %14.5e %14.5e\n',...
' ', 'Scilab',[df_true',df_scilab',err_abs_sci',err_rel_sci']);
mprintf('%s\n',strcat(repmat('-',1,120)));
