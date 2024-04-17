# lecture 2 code
using PyPlot
using BenchmarkTools

# back substitution

function back_sub(A, b)
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

function linear_solver(A, b)
    n = length(b)

    for row = 1:n
        for i = row + 1 : n
            coef = A[i, row] / A[row, row]
            for j = row : n
                A[i, j] = A[i, j] - coef * A[row, j]
            end
            # 课上的程序缺少了这行
            b[i] = b[i] - coef * b[row]
        end
    end
    return back_sub(A,b)
end

n = 200

A = zeros(n, n)
for i = 1:n
    for j = 1:n
        A[i, j] = randn(1)[1]
    end
end
b = randn(n);

A \ b
linear_solver(A, b)

@btime linear_solver(A, b)