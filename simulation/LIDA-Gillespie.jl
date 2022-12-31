using BioSimulator, Plots

NA = 6.023e23;

model = Network("LIDA")         
model <= Species("Tl", 10)      #1 Tl
model <= Species("R1", 1000)    #2 R1
model <= Species("R2", 1000)    #3 R2
model <= Species("TR1", 0)      #4 TR1
model <= Species("TR2", 0)      #5 TR2
model <= Species("Nl", 0)       #6 Nl
model <= Species("Tr", 0)       #7 Tr
model <= Species("L1", 1000)    #8 L1
model <= Species("L2", 1000)    #9 L2
model <= Species("TL1", 0)      #10 TL1
model <= Species("TL2", 0)      #11 TL2
model <= Species("Nr", 0)       #12 Nr
model <= Species("D", 0)        #13 D
model <= Species("O1", 0)       #14 O1
model <= Species("O2", 0)       #15 O2

twod = 2e10/NA;
threed = 2e7/NA;

kf_TR1   = twod;    
kf_TR2   = threed;     
kf_TL1   = twod;    
kf_TL2   = threed;     

kb_TR1   = 100;  
kb_TR2   = 100; 
kb_TL1   = 100;  
kb_TL2   = 100;   

kf_lig = 2e-2; 

kf_duplex    = twod;   
kb_duplex    = 0.37; 

kf_O1  = twod;    
kf_O2  = threed;      
kb_O1  = 100;  
kb_O2  = 100;      

model <= Reaction("fwd_1", kf_TR1, "Tl + R1 --> TR1")
model <= Reaction("bwd_1", kb_TR1, "TR1 --> Tl + R1")

model <= Reaction("fwd_2", kf_TR2, "Tl + R2 --> TR2")
model <= Reaction("bwd_2", kb_TR2, "TR2 --> Tl + R2")

model <= Reaction("fwd_3", kf_TR2, "TR1 + R2 --> Nl")
model <= Reaction("bwd_3", kb_TR2, "Nl --> TR1 + R2")

model <= Reaction("fwd_4", kf_TR1, "TR2 + R1 --> Nl")
model <= Reaction("bwd_4", kb_TR1, "Nl --> TR2 + R1")

model <= Reaction("lig_5", kf_lig, "Nl --> D")

model <= Reaction("fwd_6", kb_duplex, "D --> Tl + Tr")
model <= Reaction("bwd_6", kf_duplex, "Tl + Tr -->D")

model <= Reaction("fwd_7", kf_TL1, "Tr + L1 --> TL1")
model <= Reaction("bwd_7", kb_TL1, "TL1 --> Tr + L1")

model <= Reaction("fwd_8", kf_TL2, "Tr + L2 --> TL2")
model <= Reaction("bwd_8", kb_TL2, "TL2 --> Tr + L2")

model <= Reaction("fwd_9", kf_TL2, "TL1 + L2 --> Nr")
model <= Reaction("bwd_9", kb_TL2, "Nr --> TL1 + L2")

model <= Reaction("fwd_10", kf_TL1, "TL2 + L1 --> Nr")
model <= Reaction("bwd_10", kb_TL1, "Nr --> TL2 + L1")

model <= Reaction("lig", kf_lig, "Nr --> D")

model <= Reaction("fwd_12", kf_O1, "R1 + L1 --> O1")
model <= Reaction("bwd_12", kb_O1, "O1 --> R1 + L1")

model <= Reaction("fwd_13", kf_O2, "R2 + L2 --> O2")
model <= Reaction("bwd_13", kb_O2, "O2 --> R2 + L2")


Labels = ["Tl" "R1" "R2" "TR1" "TR2" "Nl" "Tr" "L1" "L2" "TL1" "TL2" "Nr" "D" "R" "L"]

# for time in 10 .^ (range(4,17,length=100))
#     global result = simulate(model, Direct(); tfinal = time)
# # for variable in collect(1:15)
# #     plottinio = plot(result[variable,:], summary = :trajectory,
# #         xlabel = "time", ylabel = "count", label = Labels[variable])
# #     display(plottinio)
# # end
#     display(plot(result, summary = :trajectory,
#     xlabel = "time", ylabel = "count", label = Labels))
# end

tempo = 2e10
resultout = simulate(model, Direct(); tfinal = tempo)
display(plot(resultout, summary = :trajectory, xlabel = "time", ylabel = "count", label = Labels))          


# for i = 1:length(Labels)
# plotz = plot(resultout.t,resultout[i,:],label = Labels[i])
# display(plotz)
# end
num = 13;
# display(plot(resultout.t,resultout[num,:],label = Labels[num]))


