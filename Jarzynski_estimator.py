#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 16:18:04 2023

@author: hallegot
"""
## Jarzynksi equality based free energy estimator for SMD simulations ##
## Need to have a file with 3 columns : { t, q(t), f(t) } in this order.
## we call w the work, q the position (or the CV), f the biasing force
## and F the free energy.
## TODO : add computation of the errorbars.

import numpy as np
import matplotlib.pyplot as plt
import math as m


rep = [
       "/path/to/simu1",
       "/path/to/simu2/"
      ]

filename =  [
             "data1",
             "data2"
            ]


####### ++++++ ------- Parameters ------- ++++++ #######

system = 0 #id of the system to analyze in rep and filename
n_bins_wq = 50


####### ++++++ ------- Loading and organizing files ------- ++++++ #######

time = np.loadtxt( rep[system]+filename[system], comments="#", usecols=0 )
q_t = np.loadtxt( rep[system]+filename[system], comments="#", usecols=1 )
f_t = np.loadtxt( rep[system]+filename[system], comments="#", usecols=2 )


####### ------- Separation of trajectories ------- #######
split_time = []
split_q_t = []
split_f_t = []
ntraj=0
jmin=0
dt_ref = time[1] - time[0]
dt_ref+=(dt_ref*0.01)
for j in range(len(time)-1) :
    if abs(time[j] - time[j+1]) > dt_ref :
        ntraj+=1
        split_time.append( time[jmin:j] )
        split_q_t.append( q_t[jmin:j] )
        split_f_t.append( f_t[jmin:j] )
        jmin=j+1

## Sorting function converting an array ordered by time
## in an array order by position.
## inputs : ft = time order array; x_t = time order positions;
def convert_Wt_to_Wq(ft, x_t, nptgridx) :
    ntimes = np.size( ft )
    fq = np.zeros( nptgridx )
    fq2 = np.zeros( nptgridx )
    nf = np.zeros (nptgridx, int)
    xmin = np.min(x_t)
    xmax = np.max(x_t)
    dx = (xmax-xmin)/nptgridx
    for it in range(ntimes) :
        ix = m.floor( (x_t[it] - xmin)/dx )
        if 0 < ix < nptgridx :
            fq[ix] += ft[it]
            fq2[ix] += ft[it]
            nf[ix] += 1
    fq2*=fq2
    for ix in range(nptgridx) :
        if nf[ix] == 0 :
            fq[ix] = 0.
            fq2[ix]=0.
        else :
            fq[ix] /= nf[ix]
            fq2[ix] /= nf[ix]
    std = fq2 - fq*fq
    return nf, fq, std


####### ++++++ ------- Estimator ------- ++++++ #######

## Computation of the work w(t) for each trajectory
Ntraj=len(split_q_t)
moy_expw = []
w = [[0] for _ in range(Ntraj)]
exp_w = [[1] for _ in range(Ntraj)]
start = 2
for i in range(Ntraj) :
    w_acc=0
    for j in range(start,len(split_q_t[i])-1) :
        dq = split_q_t[i][j+1] - split_q_t[i][j]
        w_acc+= 0.5*(split_f_t[i][j+1] + split_f_t[i][j]) * dq
        w[i].append( w_acc ) ## work produced by f_t(t_j) in [q(t_j+1)-q(t_j)].


## Convert w(q) into w(t)
w_q=[]
n_q=[]
w_q_std=[]
for i in range(Ntraj) :
    tmp_nf, tmp_fq, tmp_fqstd =  convert_Wt_to_Wq(w[i], split_q_t[i][start:], n_bins_wq)
    w_q.append(tmp_fq)
    n_q.append(tmp_nf)
    w_q_std.append(tmp_fqstd)

## Averaging over all the trajectories with exponentials.
w_q_array = np.array(w_q)
w_q_std = np.array(w_q_std)
w_moy = np.mean( np.exp(-1.0*w_q_array), axis=0)
w_ecty = np.std( w_q_std, axis=0)

## Computing the free energy F.
F = - np.log(w_moy)
F_ecty = np.zeros(len(w_moy))
for i in range(len(w_moy)) :
    F_ecty[i] = abs( (w_moy[i]/F[i]) * w_ecty[i] )

####### ++++++ ------- plots ------- ++++++ #######

x = np.linspace(np.min(split_q_t[0]),np.max(split_q_t[0]), len(F))

plt.title("Jarzynski equality base estimator")
plt.xlabel("position")
plt.ylabel("free energy")
plt.plot(x, F, label="Jarzynski")
plt.legend()
plt.show()
plt.close()

#np.savetxt("jarzynski", np.transpose([x, F]), delimiter='  ' )

## One source on how can be implemented such estimator.
# @article{
# doi:10.1073/pnas.1635159100,
# author = {Jeff Gore  and Felix Ritort  and Carlos Bustamante },
# title = {Bias and error in estimates of equilibrium free-energy differences from nonequilibrium measurements},
# journal = {Proceedings of the National Academy of Sciences},
# volume = {100},
# number = {22},
# pages = {12564-12569},
# year = {2003},
# doi = {10.1073/pnas.1635159100},
# URL = {https://www.pnas.org/doi/abs/10.1073/pnas.1635159100},
# eprint = {https://www.pnas.org/doi/pdf/10.1073/pnas.1635159100},
# }
