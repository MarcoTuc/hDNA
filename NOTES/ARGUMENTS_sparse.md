# Argument for why we gillespie simulate with only two initial strands
If we've got for example the following initial count:
    - 1 AAAAAA
    - 2 TTTTTT
We should also count the situation:     A A A A A A
                                  T T T T T T T T T T T T
                                             ^
                                            NICK

Which is outside the scope of the model right now.
One could argue that the formation of these structure indeed 
impacts the hybridization rates and I would agree. 
This could be object of future developments on the model. 

# TITLE
# The Hy--- Simulator: Spatially dependant Stochastic Simulation of Nucleic Acids Hybridization Kinetics