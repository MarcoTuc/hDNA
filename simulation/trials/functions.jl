#################################
####### Utility Functions #######
#################################

function printwholematrix(matrix)
    Base.print_matrix(IOContext(stdout, :limit => false), matrix)
end


###################################################
## Find the key relative to duplex concentration ##
###################################################

function keyD(ordict)
    if typeof(ordict) == OrderedDict{Symbol, Species}
        for i in range(1,size(ordict.keys)[1])
            if ordict.keys[i] == :D
                return i 
            end
        end
    elseif typeof(ordict) == Vector{Symbol}
        for i in range(1,size(ordict)[1])
            if ordict[i] == :D
                return i 
            end
        end
    end
end
    

#####################################################
## get the first passage time for duplex formation ##
#####################################################

function condition(y)
    return 1 == y
end

function get_fpt(r,t,index,duplexindex)
    index_r = findfirst(condition,r[index][:,duplexindex])
    firstpassage = t[index][index_r]
    return firstpassage
end

function fpts(r,t,d)
    # r = simulation results
    # t = vector of time vectors
    # d = duplex index 
    fpts = [get_fpt(r,t,s,d) for s in 1:size(r)[1]]
    return fpts
end


"custom functions have been included"