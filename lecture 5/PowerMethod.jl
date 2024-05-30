using LinearAlgebra

# Power Method 函数实现
function power_method(A, max_iter=1000, tol=1e-6)
    n = size(A, 1)
    q = rand(n)  # 随机初始化特征向量
    q ./= norm(q)  # 归一化
    lambda = real(dot(q, A*q))  # 初始化特征值
    
    for k in 1:max_iter
        q0 = copy(q)
        q = A * q0
        q ./= norm(q)  # 归一化
        lambda_new = real(dot(q, A*q))
        if abs(lambda - lambda_new) < tol
            break
        end
        lambda = lambda_new
    end
    
    return lambda, q
end

# 创建一个测试矩阵 A
A = [2. 3. 4. 5. 6.;
     4. 4. 5. 6. 7.;
     0. 3. 6. 7. 8.;
     0. 0. 2. 8. 9.;
     0. 0. 0. 1. 10.]

# 使用 Julia 自带的 eigen 函数计算最大模特征值和特征向量
eigen(A).values
eigenvalues, eigenvectors = eigen(A)
max_eigval_index = argmax(abs.(eigenvalues))
max_eigval = eigenvalues[max_eigval_index]
max_eigvec = eigenvectors[:, max_eigval_index]

# 运行 Power Method
lambda_power, q_power = power_method(A)

# 输出结果
println("最大模特征值 (Power Method): ", lambda_power)
println("最大模特征向量 (Power Method): ", q_power)
println("最大模特征值 (Eigen): ", max_eigval)
println("最大模特征向量 (Eigen): ", max_eigvec)