using BioSimulator
using Dictionaries
using DataStructures
using Statistics

include("functions.jl")

model = Network("AT-all")

initial = 1

model <= Species("SS",Int(initial))
model <= Species("L1",0)
model <= Species("L2",0)
model <= Species("R1",0)
model <= Species("R2",0)
model <= Species("D",0)


kf_dl = 10e7
nonspecific_correction = 1
kf_sliding = 5e6


Na = 6.023e23

kf_ssl1 = kf_dl
kb_ssl1 = 22.920792678077465

kf_ssl2 = kf_dl
kb_ssl2 = 0.007808225837362128

kf_ssr2 = kf_dl
kb_ssr2 = 0.0626785751957368

kf_ssr1 = kf_dl
kb_ssr1 = 183.99091641848813

kf_ssD = kf_dl
kb_ssD = 6.902831840829526e-08

kf_l1l2 = kf_sliding
kb_l1l2 = 17.033062396725658

kf_l2D = kf_sliding
kb_l2D = 0.4420230654574368

kf_r1r2 = kf_sliding
kb_r1r2 = 17.033062396725658

kf_r2D = kf_sliding
kb_r2D = 0.05506532191640359




model <= Reaction("ss-l1_f",  kf_ssl1, "SS --> L1")
model <= Reaction("ss-l1_b",  kb_ssl1, "L1 --> SS")

model <= Reaction("ss-l2_f",  kf_ssl2, "SS --> L2")
model <= Reaction("ss-l2_b",  kb_ssl2, "L2 --> SS")

model <= Reaction("ss-r2_f",  kf_ssr2, "SS --> R2")
model <= Reaction("ss-r2_b",  kb_ssr2, "R2 --> SS")

model <= Reaction("ss-r1_f",  kf_ssr1, "SS --> R1")
model <= Reaction("ss-r1_b",  kb_ssr1, "R1 --> SS")

model <= Reaction("ss-D_f",   kf_ssD, "SS --> D")
model <= Reaction("ss-D_b",   kb_ssD, "D --> SS")

model <= Reaction("l1-l2_f",  kf_l1l2, "L1 --> L2")
model <= Reaction("l1-l2_b",  kb_l1l2, "L2 --> L1")

model <= Reaction("l2-D_f",   kf_l2D, "L2 --> D")
model <= Reaction("l2-D_b",   kb_l2D, "D --> L2")

model <= Reaction("r1-r2_f",  kf_r1r2, "R1 --> R2")
model <= Reaction("r1-r2_b",  kb_r1r2, "R2 --> R1")

model <= Reaction("r2-D_f",   kf_r2D, "R2 --> D")
model <= Reaction("r2-D_b",   kb_r2D, "D --> R2")



runtime = 1e-2
Nmonte = 3000

include("functions.jl")

simresults = [simulate(model, Direct();tfinal=runtime) for i in 1:Nmonte]
results = [permutedims(hcat(simresults[i]...)) for i in 1:Nmonte]
timeresults = [simresults[i].t for i in 1:Nmonte]

Dindx = keyD(model.species_list.keys);
fpt1 = get_fpt(results,timeresults,1,Dindx);

mbareporcoddio = getfpts(results,timeresults,Dindx);

fpt2 = get_fptforfpts(results[1],timeresults[1],5);

M = mean(mbareporcoddio);

using PlotlyJS
plot(
    [histogram(x=mbareporcoddio, opacity=0.5, nbinsx=50)],
        Layout(
            yaxis_type="",
            xaxis_title_text="first passage time",
            yaxis_title_text="# of runs in the bin")  # put "log" for logarithmic yaxis
)

1/M