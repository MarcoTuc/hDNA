import juliacall


jl = juliacall.newmodule('Dude')

jl.seval('using BioSimulator')


model = jl.Network("AT-all")  

initial = 100

model <= jl.Species("SS",int(initial))
model <= jl.Species("L1",0)
model <= jl.Species("L2",0)
model <= jl.Species("R1",0)
model <= jl.Species("R2",0)
model <= jl.Species("D",0)


kf_dl = 4e8
nonspecific_correction = 1/5

sliding_correction = 1/2
kf_sliding = 1e7 * sliding_correction



Na = 6.023e23

kf_ssl1 = kf_dl*nonspecific_correction
kb_ssl1 = 22.920792678077465 

kf_ssl2 = kf_dl*nonspecific_correction
kb_ssl2 = 0.007808225837362128

kf_ssr2 = kf_dl*nonspecific_correction
kb_ssr2 = 0.0626785751957368

kf_ssr1 = kf_dl*nonspecific_correction
kb_ssr1 = 183.99091641848813

kf_ssD = kf_dl*nonspecific_correction
kb_ssD = 6.902831840829526e-08

kf_l1l2 = kf_sliding
kb_l1l2 = 3406.6124793451313*sliding_correction

kf_l2D = kf_sliding
kb_l2D = 88.40461309148736*sliding_correction

kf_r1r2 = kf_sliding
kb_r1r2 = 3406.6124793451313*sliding_correction

kf_r2D = kf_sliding
kb_r2D = 11.01306438328072*sliding_correction




model <= jl.Reaction("ss-l1_f",  kf_ssl1, "SS + SS --> L1")
model <= jl.Reaction("ss-l1_b",  kb_ssl1, "L1 --> SS + SS")
model <= jl.Reaction("ss-l2_f",  kf_ssl2, "SS + SS --> L2")
model <= jl.Reaction("ss-l2_b",  kb_ssl2, "L2 --> SS + SS")
model <= jl.Reaction("ss-r2_f",  kf_ssr2, "SS + SS --> R2")
model <= jl.Reaction("ss-r2_b",  kb_ssr2, "R2 --> SS + SS")
model <= jl.Reaction("ss-r1_f",  kf_ssr1, "SS + SS --> R1")
model <= jl.Reaction("ss-r1_b",  kb_ssr1, "R1 --> SS + SS")
model <= jl.Reaction("ss-D_f",   kf_ssD, "SS + SS --> D")
model <= jl.Reaction("ss-D_b",   kb_ssD, "D --> SS + SS")
model <= jl.Reaction("l1-l2_f",  kf_l1l2, "L1 --> L2")
model <= jl.Reaction("l1-l2_b",  kb_l1l2, "L2 --> L1")
model <= jl.Reaction("l2-D_f",   kf_l2D, "L2 --> D")
model <= jl.Reaction("l2-D_b",   kb_l2D, "D --> L2")
model <= jl.Reaction("r1-r2_f",  kf_r1r2, "R1 --> R2")
model <= jl.Reaction("r1-r2_b",  kb_r1r2, "R2 --> R1")
model <= jl.Reaction("r2-D_f",   kf_r2D, "R2 --> D")
model <= jl.Reaction("r2-D_b",   kb_r2D, "D --> R2")



runtime = 1e-2
Nmonte = 50000



simresults = [jl.simulate(model, jl.Direct(),tfinal=runtime) for i in range(Nmonte)]


"""
results = jl.seval('[permutedims(hcat(simresults[i]...)) for i in 1:50000]')
timeresults = [simresults[i].t for i in range(Nmonte)]

Dindx = jl.keyD(model.species_list.keys)
fpt1 = jl.get_fpt(results,timeresults,1,Dindx)

sao = jl.getfpts(results,timeresults,Dindx)

fpt2 = jl.get_fptforfpts(results[1],timeresults[1],5)

M = jl.mean(sao)


jl.plot(
    [jl.histogram(x=sao, opacity=0.5, nbinsx=5000)],
        jl.Layout(
            yaxis_type="",
            xaxis_title_text="first passage time",
            yaxis_title_text="# of runs in the bin")  # put "log" for logarithmic yaxis
)

1/M

class Simulator(object):
    pass 

import juliapkg
print("Current juliacall project folder is:", juliapkg.project())
"""