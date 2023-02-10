#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>>>>#### Simulation Exploration 

>>>>
    Decreasing minimum-nucleation has the effect of lowering average resulting rates for strands.
    I saw this by confronting 1 vs 3 minimum nucleation (lost data)
    Since my model is already under-computing rates for minimum nucleation of 3, I will need to 
    do a new simulation with some higher minimum nucleation value and see how better it behaves.
>>>>
    Did a minimum-nucleation = 4 simulation and results are very high, I need to try lowering down
    steric angles and maybe also lower down zipping one order of magnitude and see what happens. 


######################################################
>>>>>>>>>>>>>>>>>>>>>>>#### THEORETICAL CONSIDERATIONS 

<Non-Markovianity>
Right now the kinetic network we gillespie simulate is a markovian one. 
We are not inferring transitions probabilities from any other information than the current state. 
Actually this is not the case as the actual physical location of strands will influence the probability
of having this or that subsequent transition. For example let's think about a sliding going into zipping: 
'........(.+........).' sliding
'..........+..........' singlestranded
'(.........+.........)' nucleation
The sliding in this case is one like this: AAAAAAAAAG
                                                   TCTTTTTTTT 
Then it simplexes: 
Then it nucleates like: AAAAAAAAAG
                        TCTTTTTTTT
                        ^
It is clear that here there's no information about spatial permanence of objects 
which will eventually influence markovianity rendering the process actually dependant
on current state + state previous to that (it is just order 1 nonmarkovian)


<Multistrand>
The multistrand package doesn't have the minimum nucleation parameter.


10, 11, 26

