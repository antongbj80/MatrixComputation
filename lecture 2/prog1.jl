#运行按shift+enter
using PyPlot
using LinearAlgebra
using BenchmarkTools

#上三角随机生成
n = 5
A = zeros(n,n)

for i in 1:n
    for j = i:n
        A[i,j] = randn(1)[1]
    end
end
b = randn(n)

figure()
spy(A)
display(gcf())

function back_sub(A,b)
    x = zeros(n)
    for i in n:-1:1 
        x[i] = b[i]
        for j in i+1 :n
            x[i] -= A[i,j] *x[j]
        end
        x[i] = x[i] / A[i,i]
    end
    return x
end

