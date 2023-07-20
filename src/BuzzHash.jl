module BuzzHash

using SparseArrays, Random, Statistics

export sprand_fd, buzzhash, inverse

""" sprand_fd(KC, PN, p)
    
    Returns a KC by PN sparse 1/0 (technically Int8) matrix in which each entry is 1 with probability p. KC and PN abbreviate Kenyon Cells and Projection Neurons, respectively, to acknowledge the algorithm's biological source. Generally KC is much larger than PN (e.g., 2000 vs 50) and p is relatively small (e.g., 0.12). 
   
    (I.e., sprand_fd is essentially a wrapper for Base.sprand.)
    """
function sprand_fd(KC::T1, PN::T2, p::R)::SparseMatrixCSC{Int8,Int} where {T1<:Integer, T2<:Integer, R<:Real}
    tmp = sprand(Bool,KC,PN,p) # Bool to make it binary
    i,j = findnz(tmp)
    return sparse(i,j,Int8(1)) # Convert to Int8
end

""" sprand_fd(KC, PN, PN_per_KC)
    
    If the last argument is an integer sprand_fd returns a KC by PN sparse 1/0 (technically Int8) matrix in which each row has exactly PN_per_KC randomly placed 1's. This version is much slower.
    """
function sprand_fd(KC::T1, PN::T2, PN_per_KC::T3)::SparseMatrixCSC{Int8,Int} where {T1<:Integer, T2<:Integer, T3<:Integer}
    KC > 0 && PN > 0 || error("KC and PN must be positive.")   
    tmp = collect(1:PN)
    idx = 1:PN_per_KC
    ans = spzeros(Int8,KC,PN)
    for i in 1:KC
        shuffle!(tmp)
        ans[i,tmp[idx]].=Int8(1)
    end
    return ans
end

""" buzzhash(A, x, topN; clip=true)

    Return a sparse vector formed by applying A to x-mean(x), zeroizing all but the topN largest components of the result, and by default setting the remaining components to 1 ("clipping" them.) To retain the component values, call with `clip=false`. A is a random, sparse, binary (technically Int8) matrix implementing the expansion stage of the algorithm. 
"""
function buzzhash(A::SparseMatrixCSC{Int8,Int}, x::Vector{R}, topN::T; clip::Bool=true)::SparseVector{Float64,Int} where {R<:Real,T<:Integer}
    # center x and apply A
    tmp = A*(x .- mean(x))
    # find the permutation which sorts tmp in descending order
    perm = sortperm(tmp, rev=true)
    # return a sparse vector containing the largest topN elements of the result.
    return sparsevec(perm[1:topN], clip ? 1 : tmp[perm[1:topN]], size(A,1)) 
end

""" inverse(A)

    Returns a closure (function + data) which, when applied to Ax, returns x.

    When KC >> PN the PN columns of a sparse random matrix formed by sprand_fd are linearly independent with high probability. This implies that the map, y = Ax, is invertable--that x can be recovered from y. The inverse mapping is inv(A'A)A' where A' is the transpose of A and inv indicates inverse. The inversion will be at best approximate, of course, after the final hashing step in which all but the largest values are zeroized.
    """
function inverse(A::SparseMatrixCSC{Int8,Int})
    B = Matrix{Float64}(A)
    S = inv(B'*B)
    (x) -> S*(B'*x)
end


end # module
