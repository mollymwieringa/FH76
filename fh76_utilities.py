#############################################
## FH76 Stochastic Climate Model Utilities ##
#############################################

# This utility file defines the FH76 model, #
# the wind generation function, and any     #
# other associated functions for analyzing  #
# simple experiments related to climate     #
# variability.                              #

# Author: Molly M. Wieringa (2022)          #

#############################################
## Required Libraries                      ##
#############################################

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#############################################
## FH76 Model Definition                   ##
#############################################

def model(U, V, t=60*60*6, h=25, K= 0.25, feedback=False, plot=True):
    """
    Frankignoul and Hasselman's 1976 simple 
    stochastic model of climate variability.
    
    This function integrates stochastic
    atmopsheric forcing on some slab ocean 
    of prescribed mixed layer depth to 
    determine the evolution of SST anomalies
    resulting from that forcing. 

    Inputs:
    (1) U: Zonal wind
    (2) V: Meridional wind
    (3) t: timestep (in hours); default is 
           6 hours
    (4) h: ocean mixed layer depth; default
           is 25 meters
    (5) K: constant of proportionality
           between meridional wind and air-
           sea temperature difference (an
           estimate of heat transport); 
           default is 0.25 (degC/m/s)
        
    Outputs:
    (1) T: the timeseries evolution of SSTs
           in response to the forcing.
    (2) heat_flux: the actual forcing induced
           by U, V
    (3) change: the discrete time changes in 
           SST induced by the forcing.

    Author: Molly M. Wieringa (2022, Univ. of
            Washington, mmw906@uw.edu)
    """
    # Heat flux parameterizations
    C_h = 1e-3      # bulk transfer coefficient, sensible heat flux
    Bow = 3         # ratio of latent to sensible heat flux
    p_a = 1.25e-3   # density of air (g/cm3)
    p_w = 1         # density of water (g/cm3)
    Cap = 0.24      # specific heat of air (cal/g/degC)
    Cwp = 0.96      # specific heat of water (cal/g/degC)
    # K = 0.25      # const. of proportionality (temp to NS wind)

    mean_U = np.mean(np.abs(U)) 
    fbf = p_a*Cap*C_h*(1+Bow)*mean_U/(p_w*Cwp*h)

    # Model operator
    T = []
    heat_flux = []
    change = []
    for i in range(0,len(U)):
        std_f = C_h * (1+Bow) * (p_a*Cap*K*V[i]*np.abs(U[i]))/(p_w*Cwp)
        heat_flux.append(std_f)
        
        if i == 0:
            del_T = heat_flux[i]/h * t 
            T.append(del_T)
        else:
            if feedback is True:
                del_T = (heat_flux[i]/h - fbf*T[i-1]) * t
                T.append(T[i-1] + del_T)
            else:
                del_T = heat_flux[i]/h * t
                T.append(T[i-1] + del_T)
        change.append(del_T)

    if plot is True:
        # Plotting the results
        fig, axes = plt.subplots(3, 1, figsize=(8,5))

        #create boxplot in each subplot
        sns.lineplot(data=heat_flux, ax=axes[0])
        sns.lineplot(data=change, ax=axes[1])
        sns.lineplot(data=T, ax=axes[2])

        axes[0].set_title('Wind Forcing', fontweight='bold')
        axes[1].set_title('SST Change/$\Delta$t', fontweight='bold')
        axes[2].set_title('Sea Surface Temperature', fontweight='bold')

        plt.tight_layout()
    else:
        pass
       # print("If you'd like to plot the results of this single experiment, please set plot=True")

    return T, heat_flux, change

#############################################
## Generation Functions                    ##
#############################################

def generate_wind(x_corr=0.35, y_corr=0.35, len=4*10*365, mean_z_wind=5):
    """
    A function to generate a semi-realistic 
    pair of timeseries for meridional and 
    zonal winds.
    
    This function generates wind timeseries
    with a prescribed auto-correlation 
    ("memory per timestep") value by drawing
    perturbations from a random normal 
    distribution and imposing memory and a 
    mean zonal wind value."

    Inputs:
    (1) x_corr: autocorrelation of zonal
            wind; default=0.35 per Sarkar et
            al. (2002)
    (2) y_corr: autocorrelation of merid-
            ional wind; defatul = 0.35 per
            Sarkar et al. (2002)
    (3) len: number of timesteps in desired
            timeseries
    (4) mean_z_wind: mean value of the zonal
            wind component. 
        
    Outputs:
    (1) U: Zonal wind 
    (2) V: Meridional wind

    Author: Molly M. Wieringa (2022, Univ. of
            Washington, mmw906@uw.edu)
    """ 
    U = []
    V = []
    U.append(np.random.randn())
    V.append(np.random.randn())

    for i in range(0,len):
        x = mean_z_wind + x_corr * U[i] + np.random.randn()
        y = y_corr * V[i] + np.random.randn()

        U.append(x)
        V.append(y)

    U = np.array(U)
    V = np.array(V)

    return U,V

def coin_flips(number):
    """
    A function to randomly flip a coin,
    return either 1 (heads) or -1 (tails),
    and track the outcome.

    Inputs:
    (1) number: the number of coin flips to
            be performed

    Outputs:
    (2) flips: a list of the results of the 
            coin flips

    Author: Molly M. Wieringa (2022, Univ. of
            Washington)
    """
    rands = np.random.rand(number)
    flips = []
    for x in rands:
        flips.append(1 if x > 0.5 else -1)
    
    return flips