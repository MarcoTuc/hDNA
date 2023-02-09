# TODO TOMORROW:
>>> Check Kawasaki, Metropolis and Entropy/Enthalpy rate methods from Schaffer Master Thesis
>>> Get a method for histogramming first passage times of each simulation ensemble 
>>> Get a method for percent completion by simulation time 
>>> Introduce boltzmann sampling of single stranded to nucleated transition 
    > either use NUPACK integrated boltzmann sampling or
    > implement it myself by taking free energy of each nucleation and distribute it in boltzmann fashion myself
>>> Implement 6.6.2 from Schaffer Master Thesis to properly calculate K_eff 

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