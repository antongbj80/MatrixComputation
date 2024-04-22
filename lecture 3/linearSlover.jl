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

#带排序的高斯消元法（列主元）求解器
function sort_gauss_elimination(A, b)
    n = length(b)
    #排序
    for k = 1 : n-1  #从1到n-1列，最后一列不需要排序，因为就一个
        row_max = k #默认第k行k列是最大的
        for m = k+1 :n  #从k+1到n行
            #找最大行的位置
            if abs(A[m,k]) > abs(A[row_max,k])
                row_max = m
            end
        end
        #判断是否需要交换
        if row_max != k #需要交换
            # 使用临时变量temp_row来存储第k行的内容
            temp_row = A[k, :]
            temp_b = b[k]
            
            # 将第max_row行的内容复制到第k行
            A[k, :] = A[row_max, :]
            b[k] = b[row_max]
            
            # 将之前第k行（现在存储在temp变量中）的内容复制到第max_row行
            A[row_max, :] = temp_row
            b[row_max] = temp_b
        end
    end
    #消元，化上三角
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
function LDL(A,b)
    
end

#三对角矩阵求解（追赶法）器
function tri_DiagSolver(A,b)
    n = size(A, 1)
    u = zeros(n)
    y = zeros(n)

    u[1] = A[1,2]/A[1,1]
    y[1] = b[1]/A[1,1]
    for i = 2:n-1
        u[i] = A[i,i+1]/(A[i,i]-u[i-1]*A[i,i-1])
        y[i] = ((b[i]-y[i-1]*A[i,i-1])/(A[i,i]-u[i-1]*A[i,i-1]))        
    end
    y[n] = ((b[n]-y[n-1]*A[n,n-1])/(A[n,n]-u[n-1]*A[n,n-1]))
    # 回代求解x
    x = zeros(n)
    x[n] = y[n]
    for i = n-1:-1:1
        x[i] = y[i] - u[i] * x[i+1]
    end
    return x
end

