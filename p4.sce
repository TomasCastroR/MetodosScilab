// Resuelve un sistema de ecuaciones triangular superior
// Sustitucion regresiva
// Recibe la matriz de coeficientes y el vector resultado
function x = superior (A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('superior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('superior - dimensiones incompatibles entre A y b');
        abort;
    end;
    n = nA
    a = [A b]
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

// Resuelve un sistema de ecuaciones triangular inferior
// Sustitucion progresiva
// Recibe la matriz de coeficientes y el vector resultado
function x = inferior (A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('inferior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('inferior - dimensiones incompatibles entre A y b');
        abort;
    end;
    n = nA
    a = [A b]
    x(1) = a(1,n+1)/a(1,1)
    for i=2:n
        sumk = 0
        for k=1:(i-1)
            sumk = sumk + a(i,k)*x(k)
        end
        x(i) = (a(i,n+1)-sumk)/a(i,i)
    end
endfunction

// Resuelve un sistema de ecuaciones diagonal
// Recibe la matriz de coeficientes y el vector resultado 
function x = diagonal (A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('diagonal - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('diagonal - dimensiones incompatibles entre A y b');
        abort;
    end;
    x = b/diag(A)
endfunction

// Recibe la matriz de coeficientes tridiagonal y el vector resultado
function x = tridiag(A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('tridiag - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('tridiag - dimensiones incompatibles entre A y b');
        abort;
    end;
    a = [A b]
    for i=1:(nA-1)
        mji = a(i+1,i)/a(i,i)
        a(i+1,i) = 0
        a(i+1,i+1) = a(i+1,i+1)-mji*a(i,i+1)
        a(i+1,n+1) = a(i+1,n+1)-mji*a(i,n+1)
    end
    //Sustitucion regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        x(i) = a(i,n+1)- a(i,i+1)*x(i+1)/a(i,i)
    end
    
endfunction

// Eliminacion de Gauss sin pivoteo
function [x,a] = gausselim(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  
    
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A b]; // Matriz aumentada
    count = 0;
    
    // Eliminación progresiva
    n = nA;
    for i = 1:(n-1)
        for j = (i+1):n
            mjk = a(j,i)/a(i,i)
            a(j,i)=0
            a(j,(i+1):(n+mb)) = a(j,(i+1):(n+mb)) - mjk*a(i,(i+1):(n+mb))
        end
    end
    
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

// Eliminacion de Gauss sin Pivoteo
// Resuelve multiples sistemas de ecuaciones
function X = gauss_mult(A,B)
    [nA,mA] = size(A) 
    [nB,mB] = size(B)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nB then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A B]; // Matriz aumentada
    
    // Eliminación progresiva
    n = nA;
    for i = 1:(n-1)
        for j = (i+1):n
            mjk = a(j,i)/a(i,i)
            a(j,i)=0
            a(j,(i+1):(n+mB)) = a(j,(i+1):(n+mB)) - mjk*a(i,(i+1):(n+mB))
        end
    end
    X(n,1:mB) = a(n,(n+1):(n+mB))/a(n,n)
    for i = (nA-1):-1:1 do
      X(i,1:mB) = (a(i,(mA+1):(mA+mB)) - (a(i,(i+1):mA)*X((i+1):mA,1:mB)))/a(i,i)
   end 
endfunction

// Calcula la inversa de matriz. Metodo de Gauss Jordan
function X = inversa(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('inversa - La matriz A debe ser cuadrada');
        abort;
    end
    I = eye(A)
    X = gauss_mult(A,I)
endfunction

// Calculo de determinante utilizando eliminacion gaussiana
function d = determinante(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('determinante - La matriz A debe ser cuadrada');
        abort;
    end
    a = [A]
    n = nA;
    for i = 1:(n-1)
        for j = (i+1):n
            mjk = a(j,i)/a(i,i)
            a(j,i)=0
            a(j,(i+1):n) = a(j,(i+1):n) - mjk*a(i,(i+1):n)
        end
    end
    d = prod(diag(a))
endfunction

function [x,a] = gausselimPP(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana con pivoteo parcial.
    
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselimPP - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselimPP - dimensiones incompatibles entre A y b');
        abort;
    end;
    
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
        temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
        
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
    end;
    
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = a(i,i+1:n)*x(i+1:n,1)
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

// Factorizacion LU con pivoteo
// Resultado PA = LU
function [L,U,P] = factLU(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('factLU - La matriz A debe ser cuadrada');
        abort;
    end
    P = eye(A)
    L = eye(A)
    U = [A]
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:mA); U(kpivot,k:mA) = U(k,k:mA); U(k,k:mA) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        
        for i=k+1:n
            L(i,k) = U(i,k)/U(k,k)
            U(i,k:mA) = U(i,k:mA)-L(i,k)*U(k,k:mA)
        end;
    end;
endfunction

// Primero consigue la factorizacion LU con pivoteo
// Luego resuelve el sistema de ecuaciones triangulares
function [L,U,P, x] = gaussLU(A,b)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('gaussLU - La matriz A debe ser cuadrada');
        abort;
    end
    P = eye(A)
    L = eye(A)
    U = [A]
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:mA); U(kpivot,k:mA) = U(k,k:mA); U(k,k:mA) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        
        for i=k+1:n
            L(i,k) = U(i,k)/U(k,k)
            U(i,k:mA) = U(i,k:mA)-L(i,k)*U(k,k:mA)
        end;
    end;
    c = P*b
    y = inferior(L, c)
    x = superior(U, y)
endfunction

// Metodo Doolittle de factorizacion LU
function [L,U]= doolittle(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('doolittle - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(i,k)*U(k,j)
            end
            U(i,j)=A(i,j) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(m,k)*U(k,i)
                end
                L(m,i) = (A(m,i)-suma)/U(i,i)
            end
        end
    end
endfunction

// Resuelve un sistema de ecuaciones aplicando Doolittle
function [L,U,x]= d_solver(A, b)
    [L,U] = doolittle(A)
    y = inferior(L,b)
    x = superior(U,y)
endfunction
// FALTA METODO DE CROUT

function [L,U] = crout(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('crout - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(j,k)*U(k,i)
            end
            L(j,i)=A(j,i) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(i,k)*U(k,m)
                end
                U(i,m) = (A(i,m)-suma)/L(i,i)
            end
        end
    end
endfunction

function [L,U,x] = crout_solver(A,b)
    [L,U]=crout(A)
    y = inferior(L,b)
    x = superior(U,y)
endfunction

A = [1 2 3 4;1 4 9 16;1 8 27 64; 1 16 81 256]
b = [2 ;10; 44; 190]

[L,U,x] = crout_solver(A,b)
disp(A)
disp(b)
disp(L)
disp(U)
disp(L*U)
disp(A*x)

printf("doolittle\n")
[L,U,x] = d_solver(A,b)
disp(A)
disp(b)
disp(L)
disp(U)
disp(L*U)
disp(A*x)
// Factorizacion de Cholesky
function [U, ind] = CholeskyV1(A)
    eps = 1.0e-8
    n = size(A,1)
    U = zeros(n,n)
    for k = 1:n
        if k==1 then
                t = A(k,k)
        else 
                t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        end
    
        if t <= eps
            printf("Matriz no definida positiva.\n")
            ind = 0
            return
        end
        U(k,k)= sqrt(t)
        for j = k+1:n
            if k==1 then 
                        U(k,j) = A(k,j)/U(k,k)
            else 
                        U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
            end
        end
    end
    ind = 1
endfunction

// Factorizacion de Cholesky
function [U,ind] = choleskyV2(A)
    // Factorización de Cholesky.
    // Trabaja únicamente con la parte triangular superior.
    //
    // ind = 1  si se obtuvo la factorización de Cholesky.
    //     = 0  si A no es definida positiva
    //
    //******************
    eps = 1.0e-8
    //******************
    
    n = size(A,1)
    U = zeros(n,n)
    
    t = A(1,1)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(1,1) = sqrt(t)
    for j = 2:n
        U(1,j) = A(1,j)/U(1,1)
    end
        
    for k = 2:n
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        if t <= eps then
            printf('Matriz no definida positiva.\n')
            ind = 0
            return
        end
        U(k,k) = sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
        end
    end
    ind = 1

endfunction
