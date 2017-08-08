Reference: Cai, Yongyang, Kenneth L. Judd, and Jevgenijs Steinbuks (2017). A nonlinear certainty equivalent approximation method for dynamic stochastic problems. Quantitative Economics, 8(1), 117-147.

There are four GAMS code files in this package, for the examples in the following paper 
Cai, Y., K.L. Judd, and J. Steinbuks (2016). A nonlinear certainty equivalent
approximation method for dynamic stochastic problems. Quantitative Economics (forthcoming).

The code file growth4D_NLCEQ.gms is for the 2-country real business cycle problem 
in Section 4 of the above paper; Food_Energy_NLCEQ.gms solves the dynamic stochastic 
model of food and clean energy in Section 5; The code newKeynesian_NLCEQ.gms solves 
the New Keynesian DSGE model with zero lower bound in Section 6; and the code file 
growth_bind_NLCEQ.gms is for solving the RBC model with a constraint on investment 
in the supplementary. 

Running these code files require a full version of GAMS. If you do not have it, you can 
submit the code files to the online NEOS server (using the GAMS version of CONOPT): 
https://neos-server.org/neos/ for a free run. If you want to run it on a local 
demo/student version of GAMS, then you could simply reduce the time horizon used for 
NLCEQ (of course, it may reduce accuracy of solutions) to meet the limit of the 
demo/student version. For example, in growth4D_NLCEQ.gms, you can change the following line
set t time /1*51/;
to
set t time /1*11/;

 
