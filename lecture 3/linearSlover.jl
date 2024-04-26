#MyLinearSolver
#author:Tong

#下三角矩阵求解器
function forward_sub(A,b)
    x = zeros(n)
    for i = 1:n
        x[i] = b[i]
        for j = 1 : i-1
            x[i] -= A[i, j] * x[j]
        end
        x[i] = x[i] / A[i, i]
    end
    return x
end

#上三角矩阵求解器
function back_sub(A,b)
    x = zeros(n)
    for i = n:-1:1
        x[i] = b[i]
        for j = i + 1 : n
            x[i] -= A[i, j] * x[j]
        end
        x[i] = x[i] / A[i, i]
    end
    return x
end

#线性方程组求解器
function linear_solver(A,b)
    n = length(b)
    for row = 1:n
        for i = row + 1 : n
            coef = A[i, row] / A[row, row]
            for j = row : n
                A[i, j] = A[i, j] - coef * A[row, j]
            end
            b[i] = b[i] - coef * b[row]
        end
    end
    return back_sub(A,b)
end

#LU分解求解器
function LU_factorization(A,b)
    n = size(A, 1)
    L = zeros(n, n)
    U = copy(A)

    for k = 1:n
        for i = k + 1 : n
            L[i,k] = U[i,k]/U[k,k]
            for j = k : n
                U[i, j] -= L[i,k] * U[k,j]
            end
        end
    end
    for i = 1:n
        L[i,i] = 1
    end
    y = forward_sub(L,b)
    x = back_sub(U, y)
    return x
end

#LDL分解求解器  
function LDL_factorization(A,b)   
    n = size(A,1)
    L = copy(A)
    D = zeros(n,n)
    Lt = zeros(n,n)
    y = copy(b)
    x = zeros(n)

    #给L和D赋值
    for i = 1:n
        for j = i:n
            if i == j
                L[i,j] = 1
                D[i,j] = A[i,i]
            else
                L[i,j] = 0
            end
        end
    end

    #求解L和D
    for i = 1:n
        for j = 1:i-1
            for k =1:j-1
                L[i,j] = L[i,j] - D[k,k]*L[i,k]*L[j,k]
            end 
            L[i,j] = L[i,j]/D[j,j]
        end

        for j = 1:i-1
            D[i,i] -= D[j,j] * L[i,j]*L[i,j]
        end 
    end

    y = forward_sub(L,b)

    #L转置Lt
    for i = 1:n
        for j = i:n 
            if i == j
                Lt[i,j] = 1
            else
                Lt[i,j] = L[j,i]
            end
        end
    end
    #M = D-1*y
    M = zeros(n)
    for i =  1:n
        M[i] = y[i]/D[i,i]
    end
    x = back_sub(Lt,M)
    return x
end


#三对角矩阵求解（追赶法）器
function tri_DiagSolver(A,b)
    n = size(A, 1)
    u = zeros(n)
    y = zeros(n)

    u[1] = A[1,2]/A[1,1]    #1
    y[1] = b[1]/A[1,1]      #1
    for i = 2:n-1 
        u[i] = A[i,i+1]/(A[i,i]-u[i-1]*A[i,i-1])      #n-2
        y[i] = ((b[i]-y[i-1]*A[i,i-1])/(A[i,i]-u[i-1]*A[i,i-1]))  #n-2        
    end
    y[n] = ((b[n]-y[n-1]*A[n,n-1])/(A[n,n]-u[n-1]*A[n,n-1]))       #1
    # 回代求解x
    x = zeros(n)
    x[n] = y[n]    #1
    for m = n-1:-1:1           
        x[m] = y[m] - u[m] * x[m+1]          #n-1
    end
    return x
end
