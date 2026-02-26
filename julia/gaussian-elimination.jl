module Matrix
    # performs forward elimination process, forming a matrix in row echelon form,
    # function then fed into back substitution to produce a solution to a system of linear equations
    # based on pseudocode from https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode
    function gaussianElimination(matrix::Matrix{Float64})
        dimRow = size(matrix, 1)
        dimCol = size(matrix, 2)
        h = 1  # pivot row
        k = 1  # pivot column
        while h <= dimRow && k <= dimCol
            iMax = argmax(abs.(matrix[h:dimRow, k])) + h - 1
            if matrix[iMax, k] == 0
                k += 1  # no pivot in this column, index forward
            else
                # swap rows
                matrix[[h, iMax], :] = matrix[[iMax, h], :]
                for i = h+1:dimRow
                    f = matrix[i, k] / matrix[h, k]
                    # vectorized row operation (replaces inner j loop)
                    matrix[i, k+1:dimCol] = matrix[i, k+1:dimCol] - matrix[h, k+1:dimCol] * f
                    matrix[i, k] = 0  # explicitly zero out below pivot (was outside loop — bug fix)
                end
                h += 1
                k += 1
            end
        end
        # back substitution
        # matrix is assumed to be an augmented matrix [A|b] of size n x (n+1)
        n = dimRow
        x = zeros(Float64, n)
        for i = n:-1:1
            x[i] = matrix[i, n+1]  # start with the RHS value
            for j = i+1:n
                x[i] -= matrix[i, j] * x[j]
            end
            x[i] /= matrix[i, i]
        end
        return x
    end
end
