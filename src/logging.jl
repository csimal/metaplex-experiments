using Dates

function logmessage(n, error)
    time = Dates.format(now(UTC), dateformat"yyyy-mm-dd HH:MM:SS")

    maxrss = "$(round(Sys.maxrss()/1048576, digits=2)) MiB"

    logdata = (; 
        n, # iteration n
        error, # some super important progress update
        maxrss) # lastly the amount of memory being used

    println(savename(time, logdata; connector=" | ", equals=" = ", sort=false, digits=2))
end