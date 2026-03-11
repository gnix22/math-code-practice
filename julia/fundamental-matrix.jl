using LinearAlgebra
# classify absorbing states and transients
function classify_states(P)
    n = size(P, 1)
    absorbing = [i for i in 1:n if P[i, i] ≈ 1.0]
    transient = setdiff(1:n, absorbing)
    return transient, absorbing
end
# using the fundamental matrix for an absorbing Markov chain formula,
# return it as a named tuple
function fundamental_matrix(P)
    transient, absorbing = classify_states(P)
    n = length(transient)
    P_T = P[transient, transient]
    R = P[transient, absorbing]
    S = inv(I - P_T)
    return (; S, transient, absorbing)
end

function matrix_probabilities(result)
    (; S, transient) = result
    n = size(S, 1)
    F = zeros(n, n)
    for i in 1:n
        for j in 1:n
            delta    = i == j ? 1 : 0
            F[i, j]  = (S[i, j] - delta) / S[j, j]
        end
    end
    return F
end

function print_results(result, probabilities)
    (; S, transient) = result
    F = probabilities
    println("Fundamental matrix S:")
    display(round.(S, digits=4))
    println("Probabilities:")
    display(round.(F,digits=4))
end

# allows for cmdline functionality, use convention row 1/ row 2 / row 3, separated by spaces
function parse_matrix(s::String)
    rows = split(strip(s), "/")
    data = [parse.(Float64, split(strip(r))) for r in rows]
    return reduce(vcat, [permutedims(row) for row in data])
end
# default matrix for testing
P = isempty(ARGS) ? [
    0.0  0.5  0.5  0.0  0.0;
    0.2  0.0  0.3  0.5  0.0;
    0.1  0.1  0.0  0.4  0.4;
    0.0  0.0  0.0  1.0  0.0;
    0.0  0.0  0.0  0.0  1.0;
] : parse_matrix(ARGS[1])

result = fundamental_matrix(P)
probabilities = matrix_probabilities(result)
print_results(result,probabilities)
