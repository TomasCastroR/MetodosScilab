function y = Lk(x,k)
    [nx,mx] = size(x)
    r = [x(1:k-1) x(k+1:mx)]
    p = poly(r,"x","roots")
    pk = horner(p,x(k))
    y = p / pk
endfunction

function p = lagrange(x,y)
    [nx, mx] = size(x)
    p = 0
    for k=1:mx
        p = p + (Lk(x,k)*y(k))
    end 
endfunction

function d = difDividida(x,y)
    [nx,mx] = size(x)
    if mx == 1 then
        d = y(1)
    else
        d = (difDividida(x(2:mx),y(2:mx))- difDividida(x(1:mx-1),y(1:mx-1)))/(x(mx)-x(1))
    end
endfunction

function p = newton(x,y)
    [nx,mx] = size(x)
    p = y(1)
    for k=2:mx
        q = poly(x(1:k-1),"x","roots")
        dif = difDividida(x(1:k),y(1:k))
        p = p + dif*q
    end
    
endfunction

// Método de Diferencias Divididas de Newton
// Formula con multiplicaciones encajadas
function w = DD_Newton(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    s = poly(0,"x")
    n = length(x)
    w =  difDividida(x,y)
    for j=n-1:-1:1
        w = difDividida(x(1:j),y(1:j)) + (s-x(j))*w
    end
endfunction


// Error de interpolación
function w = err(p,x,cot)
    // Entrada: p = valor real, x = nodos de interpolación, cot = cota de |f^(n))|
    // Salida: w = error de interpolación en x = p
    n = length(x)
    w = cot/(factorial(n))
    for i=1:n do
        w = w*abs(p - x(i))
    end
endfunction

// ----------------- AJUSTE DE CURVA ----------------------------

// Minimos cuadrados

// Es necesario dividir el proceso de la matriz y el polinomio resultado
// ya que se debe checkear que la matriz por su transpuesta sea rango completo
function A = matriz_mc(x,n)
    [nx,mx] = size(x)
    A = ones(mx,1)
    for i=1:n
        A = [A (x')^i]
    end
endfunction

// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
function [p,err] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction



// Método de Eliminación Gaussiana con pivoteo parcial
function [x,a] = gausselimPP(A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    a = [A b]; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i=k+1:n
            if abs(a(i,k))>amax then
                kpivot = i; amax = a(i,k);
            end;
        end;
        temp = a(kpivot,:); a(kpivot,:) = a(k,:);
        a(k,:) = temp
        
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k)
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end
    end
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n)
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k)
        end
        x(i) = (a(i,n+1)-sumk)/a(i,i)
    end
endfunction

// Recibe un numero n
// Devuelve el polinomio de chebyshev de ese grado con sus raíces
function [T,r] = Chebyshev(n)
    t(1) = 1
    t(2) = poly([0],"x","r")
    for i = 3:n+1
        t(i) = poly([0 2], "x", "coeff")*t(i-1)-t(i-2) 
    end
    T = t(n+1)
    r = roots(T)
endfunction

// Raices de Chebyshev en cualquier intervalo
function x = NodosChebyshev(n,a,b)
    [pol,r] = Chebyshev(n)
    for i = 1 : n
        x(i) = ((b+a) + r(i) * (b - a))/2
    end
endfunction

// Funcion que calcula los nodos de raices segun el grado.
function r = Cheb(n)
    for k=0:n-1
        r(k+1) = cos(%pi/2*(1+2*k)/n)
    end
endfunction

// ------------------------  EJERCICIO 1 ----------------------
/*
x = [0 0.2 0.4 0.6]
y = [1 1.2214 1.4918 1.8221]
resultado = 1.395612425
v = 1/3
// Lineal
pol1 = lagrange(x(2:3),y(2:3))
printf("Aproximacion lineal con lagrange:")
disp(pol1)
disp(horner(pol1,v))

pol2 = DD_Newton(x(2:3),y(2:3))
printf("Aproximacion lineal con newton:")
disp(pol2)
disp(horner(pol2,v))

//Cubica
pol1 = lagrange(x,y)
printf("Aproximacion cubica con lagrange:")
disp(pol1)
disp(horner(pol1,v))
pol2 = DD_Newton(x,y)
printf("Aproximacion cubica con newton:")
disp(pol2)
disp(horner(pol2,v))
*/

//  --------------------------  EJERCICIO 4-----------------------
/*
x = [2,2.1,2.2,2.3,2.4,2.5]
y = [0.2239,0.1666,0.1104,0.0555,0.0025,-0.0484]

p = DD_Newton(x,y)
disp("El polinomio interpolante es: ")
disp(p)

w1 = horner(p,2.15)
err1 = err(2.15,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.15) es: "+string(w1))
disp("con error: "+string(err1)+" < 0.5D-06")

w2 = horner(p,2.35)
err2 = err(2.35,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.35) es: "+string(w2))
disp("con error: "+string(err2)+" < 0.5D-06")*/


// --------------------- EJERCICIO 5 -------------------------------
/*
// Despejamos los valores para x=0,1,2
// Utilizando los polinomios de interpolacion que son dato
// f(0) = 1, f(1) = 3 y f(2) = 3


// Despejamos f(3)
x = [1 2 3]

L1 = Lk(x,1)
L2 = Lk(x,2)
L3 = Lk(x,3)

c1 = horner(L1,2.5)*3
c2 = horner(L2,2.5)*3
c3 = horner(L3,2.5)

a = (3-c1-c2)/c3

x = [0 1 2 3]
y = [1 3 3 a]

p = lagrange(x,y)

res = horner(p,2.5) */


// -------------------- EJERCICIO 6 --------------------------
/*
// DATOS
DDF1 = 2 // f[-1]
DDF2 = 1 // f[-1,1]
DDF3 = -2 // f[-1,1,2]
DDF4 = 2 // f[-1,1,2,4]

x = [-1 1 2 4]

// y0 = DDF1
// DDF2 = (y1 - y0) / (x1 - x0) => DDF2*(x1-x0) + y0 = y1

y1 = DDF2*(x(2)-x(1))+DDF1

// DDF3 = f[x1,x2]-f[x0,x1] / (x2-x0)
// f[x1,x2] = y2 - y1 / x2 - x1
// DDF3*(x2-x0) + DDF2 = (y2 - y1) / (x2 - x1)
// [DDF3*(x2-x0)+DDF2]*(x2-x1) + y1

y2 = (DDF3*(x(3)-x(1)) + DDF2)*(x(3)-x(2)) + y1

function y = pol(x)
    y = DDF1 + (x+1)*(DDF2 + (x-1)*(DDF3 + (x-2)*DDF4))
endfunction
printf("Aproximacion de f(0): %.2f\n", pol(0))

err0 = err(0,x,33.6)

printf("Cota de error de aproximacion de f(0): %f\n", err0)*/
// El resultado es muy extremo, revisar


// -------------------------- EJERCICIO 7 -------------------------
/*
x = [0 0.15 0.31 0.5 0.6 0.75]
y = [1 1.004 1.31 1.117 1.223 1.422]

// Grado 1
A = matriz_mc(x,1)
[p1,err1] = MinCuad_pol(A,y)
disp(p1)
A = matriz_mc(x,2)
[p2,err2] = MinCuad_pol(A,y)
disp(p2)
A = matriz_mc(x,3)
[p3,err3] = MinCuad_pol(A,y)
disp(p3)

rango = [-0:0.0001:1]
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])
*/
// -------------------------------- EJERCICIO 8 -----------------------
/*
x = [4 4.2 4.5 4.7 5.1 5.5 5.9 6.3 6.8 7.1]
y = [102.56 113.18 130.11 142.05 167.53 195.14 224.87 256.73 299.5 326.72]

disp("ítem a)")

disp("(#) n=1.")
A = matriz_mc(x,1)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 1 es:")
[p1,err1] = MinCuad_pol(A,y)
disp(p1)

disp("(#) n=2.")
A = matriz_mc(x,2)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 2 es:")
[p2,err2] = MinCuad_pol(A,y)
disp(p2)

disp("(#) n=3.")
A = matriz_mc(x,3)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 3 es:")
[p3,err3] = MinCuad_pol(A,y)
disp(p3)


disp("(#) Analizamos los errores err = norm(Ax-y,2).")
disp("Para la aproximación lineal: "+string(norm(err1,2)))
disp("Para la aproximación cuadrática: "+string(norm(err2,2)))
disp("Para la aproximación cúbica: "+string(norm(err3,2)))
disp("Podemos decir que es mejor la aproximación cúbica en este caso.")

inf = min(x)-1
sup = max(x)+1
rango = [inf:0.001:sup]
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])
*/

// ----------------------------- EJERCICIO 9 -------------------------------
/*
function y=f(x)
    y = 1./(1+x.^2)
endfunction

sup = 5
inf = -sup
rango = [inf:0.0001:sup]
func = f(rango)

for n=2:2:14
    if (n <> 8 && n <> 12)
        nodos = linspace(inf,sup,n+1)
        p = lagrange(nodos,f(nodos))
        aprox = horner(p,rango)
        plot(rango,aprox,"b")
        plot(rango, func,"r")
        plot(rango, abs(aprox-func),"g")
        h1 = legend(["Aprox Grado "+ string(n), "1/1+x^2", "Error"])
        figure
        sleep(1,"s")
    end
end*/

// ------------------------------- EJERCICIO 10 ------------------------------
/*
rango = [-1:0.0001:1]
func = exp(rango)

nodos = Cheb(4)
p = DD_Newton(nodos',exp(nodos'))


plot(rango,horner(p,rango),"r")
plot(rango,func,"b")
plot(rango,abs(horner(p,rango)-func),"g")
h1 = legend(["Aprox Grado 3", "exp(x)", "Error"])
a=get("current_axes")//get the handle of the newly created axes
a.axes_visible="on"; // makes the axes visible
a.y_location= "middle"*/

// ---------------- EJERCICIO 11 -----------------------------
/*
a=0
b=%pi/2
rango = [a:0.00001:b]
func = cos(rango)
grado = 3

nodos = NodosChebyshev(grado+1,a,b)

p = DD_Newton(nodos',cos(nodos'))

plot(rango,horner(p,rango),"r")
plot(rango,func,"b")
plot(rango,abs(horner(p,rango)-func),"g")
h1 = legend(["Aprox Grado 3", "cos(x)", "Error"])*/




