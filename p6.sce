function [radios, centros] = cotasGers(A)
    [n,m] = size(A)
    centros = diag(A)
    radios = sum(abs(A),'c') - abs(centros)
    
    // Imprimo las cotas
    for i = 1:n
        printf("Un autovalor se puede encontrar en el círculo de centro %f y radio %f\n", centros(i), radios(i))
    end
endfunction

function A = Aeps(eps)
    A = [1 -1 0; -2 4 -2; 0 -1 1+eps]
endfunction


function circ(r,x,y)
    xarc(x-r,y+r,2*r,2*r,0,360*64)
endfunction

function [radios, centros] = Gers(A)
    [n,m] = size(A)
    centros = diag(A)
    radios = sum(abs(A),'c') - abs(centros)
    
    mx = round (min(centros - radios)) - 1
    my = round (min(-radios)) - 1
    
    Mx = round (max(centros + radios)) + 1
    My = round (max(radios)) + 1
    
    rect = [mx,my,Mx,My]
    plot2d(0,0,-1,"031","",rect)
    xgrid()
    for i=1:n
        circ(radios(i),centros(i), 0)
    end
endfunction

function [radios, centros] = CircGersValor(A)
    [n,m] = size(A)
    centros = diag(A)
    radios = sum(abs(A),'c') - abs(centros)
    
    mx = round (min(centros - radios)) - 1
    my = round (min(-radios)) - 1
    
    Mx = round (max(centros + radios)) + 1
    My = round (max(radios)) + 1
    
    rect = [mx,my,Mx,My]
    eigen = spec(A)
    plot2d(real(eigen),imag(eigen),-3,"031","",rect)
    xgrid(2,1,7)
    for i=1:n
        circ(radios(i),centros(i), 0)
    end
endfunction

function [v, zn, iters]= mpotencia(A,z0,eps,maxiter)
    V = max(abs(spec(A)))
    v = 0
    iters = 1
    w = A*z0
    zn = w/norm(w)
    [m,j] = max(abs(w))
    v = m/z0(j) // Falta checkear que z0(j) != 0
    er1 = abs(V-v)
    er2 = norm(zn-z0)
    err = max(er1,er2)
    while (iters <= maxiter && err > eps)
        z0 = zn
        w = A*z0
        zn = w/norm(w,'inf')
        [m,j] = max(abs(w))
        v = m/z0(j) // Falta checkear que z0(j) != 0
        er1 = abs(V-v)
        er2 = norm(zn-z0,'inf')
        err = max(er1,er2)
        iters = iters+1
    end
    printf("iteraciones: %d\n", iters)
endfunction
    
    
// --------------------- EJERCICIO 1 ----------------------------

A = [1 0 0; -1 0 1; -1 -1 2]
B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
D = [4 -1 0; -1 4 1; -1 -1 4]
E = [3 2 1; 2 3 0; 1 0 3]
F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]

// A
printf ("Matriz A\n\n")
cotasGers(A)
cotasS = spec(A)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)

// B
printf ("Matriz B\n\n")
cotasGers(B)
cotasS = spec(B)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)

// C
printf ("Matriz C\n\n")
cotasGers(C)
cotasS = spec(C)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)

// D
printf ("Matriz D\n\n")
cotasGers(D)
cotasS = spec(D)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)
// E
printf ("Matriz A\n\n")
cotasGers(E)
cotasS = spec(E)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)

// F
printf ("Matriz F\n\n")
cotasGers(F)
cotasS = spec(F)
printf("Autovalores calculos con la funcion spec:")
disp(cotasS)
    
    
//  --------------------- EJERCICIO 3 --------------------------
eps = 0.1

for i=1:10
    printf("Epsilon = %f\n", eps*i)
    A = Aeps(eps*i)
    p = poly(A,"x")
    printf("Polinomio caracteristico:")
    disp(p)
    printf("Raices del polinomio caracteristico:")
    disp(roots(p))
    printf("Autovalores de A:")
    disp(spec(A))
    printf("\n")
end
// ----------------- EJERCICIO 4 --------------------------

//A = [4 1 0; 1 4 1; 0 1 4]
//Gers(A)
//CircGersValor(A)

// ----------------- EJERCICIO 5 -----------------
/*
A1 = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
A2 = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
Z0 = ones(4,1)
eps = 10^(-8)
maxiter = 1000
[v1,z1, iters1] = mpotencia(A1,Z0,eps,maxiter)
printf("Matriz A1\n")
printf("Autovalor dominante: %f\n", v1)
printf("Autovector asociado:")
disp(z1)
printf("\n")

[v2,z2, iters2] = mpotencia(A2,Z0,eps,maxiter)
printf("Matriz A2\n")
printf("Autovalor dominante: %f\n", v2)
printf("Autovector asociado:")
disp(z2)*/
