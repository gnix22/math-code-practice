# performs forward elimination process, forming a matrix in row echelon form,
# function then fed into back substitution to produce a solution to a system of linear equations
# based on pseudocode from https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode
function gaussianElimination(listIn, dimRow, dimCol)
    matrix = reshape(listIn, dimRow, dimCol)
    h = 1 # initialize the pivot row
    k = 1 # initialization of pivot column
    while h <= dimRow && k <= dimCol
        # find the kth pivot 
        iMax = h
        for i = h:dimRow
            if abs(matrix[i,k]) > abs(matrix[iMax,k])
                iMax = i
        end
        if matrix[iMax,k] == 0
            k += 1 # no pivot in this column, move on
        else
            #swap rows
            matrix[[h,iMax], :] = matrix[[iMax,h], :]
            for i = h+1:dimRow
                f = matrix[i, k]/matrix[h,k]
                # fill column below pivot with zeros
                matrix[i,k]=0
                for j=n+1:dimCol
                    matrix[i,j]=matrix[i,j]-matrix[h,j]*f
                end
            end
            h+=1
            j+=1
        end
    end

