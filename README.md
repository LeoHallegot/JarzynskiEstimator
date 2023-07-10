# JarzynskiEstimator
Simple free energy estimator based on Jarzynski equality for biased molecular dynamics simulations.

This code has been made during my M2 Internship in the Laboratoire de Biochimie Théorique (LBT-IBPC) at Paris.
It uses Jarzynski's equality $e^{\Delta F /k_BT} = e^{-w / k_BT}$ to estimate free energy over a set of biased MD trajectories.
The code has been implemented for a single file containing all the trajectories concatenated. It needs to have the time, position and biasing force ordered in three columns.
