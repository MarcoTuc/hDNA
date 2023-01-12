#######################################################
#### About the simulator

!) Implement the secondary structure dependance 
    - implement it in Strand class
    - implement it upwards to: 
        -- Complex 
        -- Chamber 
        -- Kinetwork 
        -- Simulator

!!) Implement a method inside the Kinetwork class 
    to help produce visualizations of the kinetwork
    (see helpernotes4 for how to implement it)
    - get some coloring stuff to divide 
        -- offnucleations - blue
        -- slidings       - dark blue 
        -- onnucleations  - green
        -- zippings       - dark green
        -- singlestranded - yellow
        -- duplexed       - red 

!!!!) Add kinetic constants to edge creation routines in Kinetwork class
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
    

!!!) Implement the Simulator from the Kinetwork 
    - Translate simhelper1.ipynb code into the class 


#######################################################
#### Experimental verification

!) Make a routine for simulating the model over hertel's data 