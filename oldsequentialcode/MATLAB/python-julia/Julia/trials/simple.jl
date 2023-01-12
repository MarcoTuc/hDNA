using BioSimulator
using Plots

model = Network('m')

initial = 20

model <= Species("X",Int(initial))
model <= Species("Y",Int(initial))
model <= Species("Z",0)

kf = 2e7
kb = 2

model <= Reaction("forward",kf,"X + Y --> Z")
model <= Reaction("backward",kb,"Z --> X + Y")

runtime = 1e-8
result = simulate(model, Direct();tfinal=runtime)

plot(result, summary = :trajectory, xlabel = "time", ylabel = "count")
e = 1;
ht = findall(x -> initial/2 - e < x < initial/2 + e, result[1,:])
vline!(result.t[ht])
halftime = result.t[ht]

