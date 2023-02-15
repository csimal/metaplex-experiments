using DelimitedFiles
using Graphs


function contiguous_usa()
    es, _ = readdlm("networks/edgesUSA.csv", ',', Int, header=true)
    es .+= 1
    return  SimpleGraphFromIterator([Edge(es[i,1] => es[i,2]) for i in 1:size(es,1)])
end

function london_transport()
    es = readdlm("networks/LondonTransportEdgesSimplified.csv", ',', Int)
    es .+= 1
    return  SimpleGraphFromIterator([Edge(es[i,1] => es[i,2]) for i in 1:size(es,1)])
end