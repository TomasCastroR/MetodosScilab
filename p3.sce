//Recibe una funcion continua en un intervalo y dos puntos del intervalo
//Convergencia lineal asegurada

function x= Biseccion(f,a,b,tol)
    if f(a)*f(b) > 0 then
        x = %nan;
    else
        c = (a+b)/2;
        while abs(b-c) > tol
            if f(b)*f(c) <= 0 then
                a = c;
            else
                b = c;
            end
         c = (a+b)/2;
     end
     x=c
    end
endfunction


//Recibe una funcion y dos puntos de un intervalo
//Convergencia mas rapida que lineal
//Convergencia no asegurada
//Tratar de que F'(a) != 0

function x = Secante(f,x0,x1,tol)
    actual = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
    anterior = x1;
    while abs(actual-anterior) > tol
        siguiente = actual - f(actual) * (actual-anterior)/(f(actual)-f(anterior))
        anterior = actual
        actual = siguiente
    end
    x = actual
endfunction

function x = PuntoFijo(g,x0, tol)
    actual = x0
    anterior = actual
    while abs(actual - anterior) > tol then
        anterior = actual
        actual = g(actual)
    end
    x = actual
endfunction


function x = metodoNewton(f, x0, eps)
    anterior = x0
    actual = x0 - f(x0)/numderivative(f, x0)
    while (abs(actual - anterior) > eps)
        anterior = actual
        actual = anterior - f(anterior)/numderivative(f, anterior)
    end
    x = actual
endfunction

//Recibe una funcion y un punto del espacio vectorial
//El punto debe ser una estimacion cercana de la raiz
//Convergencia cuadratica
//Convergencia no asegurada
//Asegurarse que F'(a) != 0 

function y = MetodoNewtonMult(f,x0,tol)
    anterior = x0
    jacobiana = numderivative(f,x0)
    y = a - (inv(jacobiana))*f(x0)
    while norm(y-anterior) > tol
        jacobiana = numderivative(f,anterior)
        anterior = y
        y = anterior - (inv(jacobiana))*f(anterior)
    end
endfunction

//Recibe una funcion y dos puntos del intervalo
//Convergencia asegurada

function x = RegulaFalsi(f,a,b,tol)
    c = b - f(b)*(b-a)/(f(b)-f(a))
    while f(c) > tol
        if(f(a)*f(c) < 0) then
            b = c
        else
            a = c
        c = b - f(b)*(b-a)/(f(b)-f(a))
        end
    end
    x = c
endfunction
