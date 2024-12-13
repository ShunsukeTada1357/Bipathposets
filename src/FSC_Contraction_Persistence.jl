<<<<<<< HEAD
#We contract birth and death of topological features [b,d] to some [b',d']. 
function contractbirth(birth,FSC)# I is an interval [b,d]
    return FSC[1][birth][2]
end

function contractdeath(death,FSC)
    d = death+1
=======
#We contract birth and death of topologival features [b,d] to some [b',d']. 
function whenbirthdaysFSC(FSC)
    n = FSC[2]
    days = []
    dd = 0
    for i in 1:n
        a = length(filter(x->x[2]==i,FSC[1]))
        dd = a + dd
        push!(days, dd)
    end
    return days
end

function contractbirth(I,days)# I is an interval [b,d]
    b = I[1]
    for i in 1:length(days)
        if b <= days[i]
            return i
        end
    end
    return length(days)
end

function contractdeath(I,FSC)
    d = I[2]+1
>>>>>>> origin/dev
    if d> length(FSC[1])
        return FSC[2]
    else
        return FSC[1][d][2]-1
    end
end

function contractinterval(I, FSC)
<<<<<<< HEAD
    return [contractbirth(I[1],FSC), contractdeath(I[2],FSC)]
end
=======
    days =whenbirthdaysFSC(FSC)
    return [contractbirth(I,days), contractdeath(I,FSC)]
end
>>>>>>> origin/dev
