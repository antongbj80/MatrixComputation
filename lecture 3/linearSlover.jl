#MyLinearSolver
#author:AnTong

#下三角矩阵求解
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

#上三角矩阵求解
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

#线性方程组求解
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

#LU分解求解
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

#LDL分解求解
function LDL(A,b)
    
end

#三对角矩阵求解
function tri_DiagSolver(A,b)
    
end

#随机生成n阶矩阵
n = 100
A = zeros(n, n)
for i = 1:n
    for j = 1:n
        A[i, j] = randn(1)[1]
    end
end
b = randn(n);
x = LU_factorization(A,b)
println(x)
A\b