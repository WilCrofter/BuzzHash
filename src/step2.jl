
function sprand_fd{T<:Integer}(KC::T, PN::T, PN_per_KC::T)
    tmp = spzeros(Bool,PN)
    for i in 1:PN_per_KC tmp[i]=true end
    ans = spzeros(Bool,KC,PN)
    for i in 1:KC ans[i,:]=shuffle(tmp) end
    return ans
end

function sprand_fd{T<:Integer, R<:Real}(KC::T, PN::T, density::R)
    sprand(Bool,KC,PN)
end
