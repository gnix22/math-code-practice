# disjoint cycle decomposition algorithm for permuation groups
# based on the proposition that states: Every permutation of a finite set can be written as a cycle or as a product of disjoint cycles.
# algorithm idea obtained from https://www.homepages.ucl.ac.uk/~ucahmto/0005_2021/Ch2.S14.html
function disjointCycleDecomposition(p::Array{Int64})
    for i=1:size(p)
        if p[i]==i
            println(i, "goes to ", p[i])
        # otherwise increment.. need to figure out the mapping logic.
