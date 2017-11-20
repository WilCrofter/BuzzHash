module BuzzHash

""" sprand_fd(KC, PN, PN_per_KC)
    
    Returns a KC by PN sparse 1/0 (technically Bool)  matrix in which each row has exactly PN_per_KC randomly placed 1's. KC and PN abbreviate Kenyon Cells and Projection Neurons, respectively, to acknowledge the algorithm's biological source.
    """
function sprand_fd{T1<:Integer, T2<:Integer, T3<:Integer}(KC::T1, PN::T2, PN_per_KC::T3)::SparseMatrixCSC{Bool,Int}
    KC > 0 && PN > 0 || error("KC and PN must be positive.")   
    tmp = spzeros(Bool,PN)
    for i in 1:PN_per_KC tmp[i]=true end
    ans = spzeros(Bool,KC,PN)
    for i in 1:KC ans[i,:]=shuffle(tmp) end
    return ans
end

""" sprand_fd(KC, PN, p)
    
    If the last argument, p, is a real number between 0.0 and 1.0, a sparse matrix whose entries are 1 (true) with probability, p. In other words, sprand_fd is an alias for Base.sprand.
    """
function sprand_fd{T1<:Integer, T2<:Integer, R<:Real}(KC::T1, PN::T2, p::R)::SparseMatrixCSC{Bool,Int}
    return sprand(Bool,KC,PN,p)
end

""" buzzhash(A, x, topN)

    Return a sparse vector formed by applying A to x-mean(x), and zeroizing all but the topN largest components of the result. A is a random, sparse, binary (technically Bool) matrix implementing the expansion stage of the algorithm. 
"""
function buzzhash{R<:Real,T<:Integer}(A::SparseMatrixCSC{Bool,Int}, x::Vector{R}, topN::T)::SparseVector{Float64,Int}
    # center x and apply A
    tmp = A*(x-mean(x))
    # find the permutation which sorts tmp in descending order
    perm = sortperm(tmp, rev=true)
    # return a sparse vector containing the largest topN elements of the result.
    return sparsevec(perm[1:topN], tmp[perm[1:topN]]) 
end

end # module
