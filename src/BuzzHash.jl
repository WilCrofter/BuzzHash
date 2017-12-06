module BuzzHash

export sprand_fd, buzzhash, clip, clip!, inverse

""" sprand_fd(KC, PN, PN_per_KC)
    
    Returns a KC by PN sparse 1/0 (technically Bool)  matrix in which each row has exactly PN_per_KC randomly placed 1's. KC and PN abbreviate Kenyon Cells and Projection Neurons, respectively, to acknowledge the algorithm's biological source. Generally KC is much larger than PN (e.g., 2000 vs 50) and PN_per_KC is relatively small (e.g., 5% of KC). 
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

""" buzzhash(A, x, topN; clip=true)

    Return a sparse vector formed by applying A to x-mean(x), zeroizing all but the topN largest components of the result, and by default setting the remaining components to 1 ("clipping" them.) To retain the component values, call with `clip=false`. A is a random, sparse, binary (technically Bool) matrix implementing the expansion stage of the algorithm. 
"""
function buzzhash{R<:Real,T<:Integer}(A::SparseMatrixCSC{Bool,Int}, x::Vector{R}, topN::T; clip::Bool=true)::SparseVector{Float64,Int}
    # center x and apply A
    tmp = A*(x-mean(x))
    # find the permutation which sorts tmp in descending order
    perm = sortperm(tmp, rev=true)
    # return a sparse vector containing the largest topN elements of the result.
    return sparsevec(perm[1:topN], clip ? 1 : tmp[perm[1:topN]], size(A,1)) 
end

""" inverse(A)

    When KC >> PN the PN columns of a sparse random matrix formed by sprand_fd are linearly independent with high probability. This implies that the map, y = Ax, is invertable--that x can be recovered from y. The inverse mapping is inv(A'A)A' where A' is the transpose of A and inv indicates inverse. The inversion will be at best approximate, of course, after the final hashing step in which all but the largest values are zeroized.
    """
function inverse(A::SparseMatrixCSC{Bool,Int})::Matrix{Float64}
    nz = findnz(A)
    B = sparse(nz[1],nz[2],ones(Float64,length(nz[3])))
    tmp = Matrix(B'*B)
    rank(tmp) >= size(A,2) || error("Bad luck. The columns of A are not independent. Try a new A.")
    return inv(tmp)*B'
end


end # module
