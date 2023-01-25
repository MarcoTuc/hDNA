#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#### MODELING 

# !!! DONE: fixed parameters 
Add kinetic constants to edge creation routines in Kinetwork class
        we have the following forward rates:
            - nucleations
                this gets divided by the total number of off and on register nucleations possible 
            - sliding
                for now this is a constant value but I want to make it 
                dependant on the "separation" between one sliding and 
                the next but I don't really know how atm. 
            - zipping
                for now this is just going to be the diffusion limited rate
                without any geometric constraints

# !!! DONE: see simulator.py 
Implement the Simulator from the Kinetwork 
    - Translate simhelper1.ipynb code into the class 

# ! DONE: useless, I'm using gephi in the end 
Implement a method inside the Kinetwork class 
    to help produce visualizations of the kinetwork
    (see helpernotes4 for how to implement it)
    - get some coloring stuff to divide 
        -- offnucleations - blue
        -- slidings       - dark blue 
        -- onnucleations  - green
        -- zippings       - dark green
        -- singlestranded - yellow
        -- duplexed       - red 


#######################################################
>>>>>>>>>>>>>>>>>>>>>>>>>#### Experimental verification

!) Make a routine for simulating the model over hertel's data 

!) Put numbers in front of gext graph names when saving for order 


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