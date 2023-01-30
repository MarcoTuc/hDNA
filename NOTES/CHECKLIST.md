# TODO TOMORROW:
>>>> 1. Run simulations with parameters: zipping 1e8, sliding 5e6, minnucleation 2 and see what happens. 
>>>> 2. Refine the sliding mechanism as detailed below 
>>>> 3. Understand why strands with no kinetic trapping slidings exhibit such high rates. 

# Slidings
Slidings happen to be kinetic traps because they are: 
1.  Very fast to go to the next configuration
2.  If next configuration is very unstable they are even faster at going back to the conf before 
3.  But also the trapping configuration will always want to go back to the other configuration
    since going back to the simplex is even less kinetically likely. I need to avoid these kinetic traps.  
Need to implement some kind of fraying mechanism expecially for these trapped zippings. 
Also there is the situation happening where some zippings don't have a fast forward zipping but just are 
stable enought that the simulation traps there (check nosliding_1 simulations)




# DONE
Fix the motherfucking inverted zipping rate constants why the fuck are they there 