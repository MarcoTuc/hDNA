#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#### INVESTIGATION
>>>>
    Some strands are systematically getting over/under evaluated by almost one order of magnitude. I need to check in detail what's happening in there, if there are any buggy things taking place that need to be fixed. 

#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#### MODELING

# ! 
Implement the secondary structure dependance 
    - implement it in Strand class
    - implement it upwards to: 
        -- Complex 
        -- Chamber 
        -- Kinetwork 
        -- Simulator


# !!!! 
Dimensional analysis of constants and formulas inside (DONE A FIRST CHECK, DO A SECOND CHECK IF NEEDED)
        Kinetics class 
            -   there's clearlys some wrong shit as one can see
                by the orders of magnitude you get for diffusion


#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>#### Experimental verification



#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>#### IMPROVEMENTS & BUGS TO FIX

# !!!
# Sliding Factor
Right now the sliding factor is computed as: 
sf = 1/(all possible nucleations)
But actually what I should do is: 
sf = 1/(all possible slidings + 1) 
where the 1 refers to the native sliding 
Then for zippings take a zipping factor
zf = 1/(all possible zippings)
>>>> Check line 94 in simulator.py


# !!!! 
# Off-register nonzipping
By looking at off-register nucleations we can see that there are direct nucleations going on directly from singlestranded to lots of base-pairs.
This entails that off-register nucleations will have a profound impact on the dynamics since there are very short pathways connecting single stranded to duplex like:

        (example is made with few nucleotides long strands as an example)
        (say that allowed minimum nucleation is 3 base pairs long and we have octamers)

        the following pathway is possible:

        ........+........
        .(((((((+.)))))))
        ((((((((+))))))))

        which is unrealistic. 

------- How to solve this 
            1\  Cut off-register nucleations with a number of base pairs longer than 
                some value plus the minimum nucleation length. Add some value because otherwise
                we're limiting the trajectory too much. 
            2\  When off-register collisions entailing the presence of more than minimum nucleation
                number of base pairs happen, model them as:
>>>>>>>>>>>>    ss + ss --> off-nucleation --> zipping --> sliding --> duplex 
        
            I think that solution 2 is the physically realistic one but it needs me to do some 
            quite annoying tweaks to all the classes Strand, Complex and Chamber
            Solution 1 on the other hand is kinda physically unrealistic but at least it is easy
            to implement since it will just mean to cut some edges out of the kinetwork graph 

------- What to do about it
            1\  I will first continue implementing this stuff as it is with the short pathway 
                and see what happens to the results
            2\  Then I will first try the solution number 1 and check experiments
            3\  If I'll have time left I will go for solution number 2 
                to make things more realistic (hopefully more accurate)

      
# !!!! Nonsense inside slidings thermodynamics!   
Basal sliding rate is 2e7.
What can happen is that the free energy difference between two slidings is positive 
such that the fws is 2e7 and the bwd is some orders of magnitude over 2e7 (say 2e9).
How do I fix this? 

# !!!! Wrong zipping sides: 
Some forward zippings are actually referring to backward zippings. 
Fix this (see wrongdoingszippingline64.html line 64)
In particular it is zipping_55 that has fwd and bwd switched
>>>> By looking at new trajectory plots I see that this bug is happening quite often. 

         




    