using LinearAlgebra

function CG(A, b, x0, eps, max_iter)
    x = copy(x0)    
    r = b - A*x   # 初始残差r
    d = copy(r)     #初始的搜索方向
    r0 = dot(r,r)  #dot点积
    if sqrt(r0) < eps  #判断初始的残差
        iter = 0
        println("The initial residuals have reached the convergence condition")
        return x, iter  # 如果满足，返回初始解和迭代次数 0
    end
    for num in 1:max_iter
        iter = num
        Ad = A*d
        alpha = r0 / dot(d,Ad)
        x = x + alpha*d
        r = r - alpha*Ad
        rdot = dot(r,r)
        println("$num of residuals = ", sqrt(rdot))
        if sqrt(rdot) < eps
            println("Residuals reaching the convergence condition")
            break
        end
        beta = rdot / r0
        d = r + beta * d
        r0 = rdot
    end
    return x , iter
end


A = [4. -1 -1 -1;
     -1 4. -1 -1;
     -1 -1 4. -1;
     -1 -1 -1 4.]
b = [10, 15, 13, 12]
A\b
x0 = zeros(size(A, 1))
eps = 1e-33
max_iter = 200
x, iter = CG(A, b, x0, eps, max_iter)

println("Solution_x: ", x)  # 打印解向量
println("Iterations_num: ", iter)  # 打印迭代次数