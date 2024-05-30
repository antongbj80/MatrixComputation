#QR分解，基于Givens旋转实现
using LinearAlgebra

function qr_algorithm_hessenberg(H)

    n = size(H, 1)
    Q = Matrix(I, n, n)  # 初始化Q为单位矩阵
    R = copy(H)           # 初始化R为H的副本

    for m in 1:n-1
        # 构造Givens旋转矩阵
        for i in m+1:n-1
            if abs(R[i-1, i-1]) < 1e-10
                continue
            end
            r = hypot(R[i-1, i-1], R[i, i-1])
            cs = R[i-1, i-1] / r
            sn = R[i, i-1] / r
            R[i-1, i-1] = r
            R[i, i-1] = 0.0

            # 应用Givens旋转到R
            for j in i-1:n
                R[i-1, j], R[i, j] = cs * R[i-1, j] - sn * R[i, j], sn * R[i-1, j] + cs * R[i, j]
            end

            # 应用Givens旋转到Q
            for k in 1:n
                Q[k, i-1], Q[k, i] = cs * Q[k, i-1] - sn * Q[k, i], sn * Q[k, i-1] + cs * Q[k, i]
            end
        end
    end

    return Q, R
end

# 示例：创建一个上Hessenberg矩阵
H = [1.0 2.0 3.0 0.0
     0.0 1.0 2.0 3.0
     0.0 0.0 1.0 2.0
     0.0 0.0 0.0 1.0]

# 调用QR算法
Q, R = qr_algorithm_hessenberg(H)

# 输出结果
println("Q:\n", Q)
println("R:\n", R)

#与自带的QR分解对比
Q, R = qr(H)
println("Q(qr):\n", Q)
println("R(qr):\n", R)