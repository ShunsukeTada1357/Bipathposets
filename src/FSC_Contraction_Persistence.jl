#We contract birth and death of topological features [b,d] to some [b',d']. 
function contractbirth(birth,FSC)# I is an interval [b,d]
    return FSC[1][birth][2]
end

function contractdeath(death,FSC)
    d = death+1

    if d> length(FSC[1])
        return FSC[2]
    else
        return FSC[1][d][2]-1
    end
end

function contractinterval(I, FSC)
    return [contractbirth(I[1],FSC), contractdeath(I[2],FSC)]
end