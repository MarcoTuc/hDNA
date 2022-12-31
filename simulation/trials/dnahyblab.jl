using BioSimulator
using Plots
using Dictionaries
using DataStructures

include("functions.jl")

model = Network("AT-all")

initial = 2

model <= Species("SS",Int(initial))
model <= Species("L3",0)
model <= Species("L5",0)
model <= Species("R5",0)
model <= Species("D",0)


kf_dl = 2e6
nonspecific_correction = 1
kf_sliding = 5e3


Na = 6.023e23

kf_ssl3 = kf_dl/Na
kb_ssl3 = 3.220377e+03/Na

kf_ssl5 = kf_dl/Na
kb_ssl5 = 1.023724e+01/Na

kf_ssr5 = kf_dl/Na
kb_ssr5 = 1.810685e+02/Na

kf_l3l5 = kf_sliding/Na
kb_l3l5 = 1.589447e1/Na

kf_l5D = kf_sliding/Na
kb_l5D = 1.040815e+00/Na

kf_r5D = kf_sliding/Na
kb_r5D = 2.034481e-03/Na

kf_ssD = kf_dl/Na
kb_ssD = 2.034481e-03/Na


model <= Reaction("ss-l3_f",  kf_ssl3, "SS + SS --> L3")
model <= Reaction("ss-l3_b",  kb_ssl3, "L3 --> SS + SS")

model <= Reaction("ss-l5_f",  kf_ssl5, "SS + SS --> L5")
model <= Reaction("ss-l5_b",  kb_ssl5, "L5 --> SS + SS")

model <= Reaction("ss-r3_f",  kf_ssr5, "SS + SS --> R5")
model <= Reaction("ss-r3_b",  kb_ssr5, "R5 --> SS + SS")

model <= Reaction("l3-l5_f",  kf_l3l5, "L3 --> L5")
model <= Reaction("l3-l5_b",  kb_l3l5, "L5 --> L3")

model <= Reaction("l5-D_f",   kf_l5D, "L5 --> D")
model <= Reaction("l5-D_b",   kb_l5D, "D --> L5")

model <= Reaction("r5-D_f",   kf_l5D, "R5 --> D")
model <= Reaction("r5-D_b",   kb_l5D, "D --> R5")

model <= Reaction("ss-D_f",   kf_ssD, "SS + SS --> D")
model <= Reaction("ss-D_b",   kb_ssD, "D --> SS + SS")


runtime = 1e25
Nmonte = 1200

simresults = [simulate(model, Direct();tfinal=runtime) for i in 1:Nmonte]
results = [permutedims(hcat(simresults[i]...)) for i in 1:Nmonte]
timeresults = [simresults[i].t for i in 1:Nmonte]

Dindx = keyD(model.species_list.keys);
fpt1 = get_fpt(results,timeresults,5,Dindx)

#TODO
# compute avg and variance of fpt1