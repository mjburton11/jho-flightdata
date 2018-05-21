" process flight data from raspberri pi"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir, mkdir, rename
from processdata import parse_pidata, plot_params, trim_data, print_flightstats

kts2ms = 0.5144444
lbs2N = 4.44822
ml2m3 = 1e-6
min2sec = 1/60.
ft2m = 0.3048
sm2kgkWhr = 1000.*60.*60.
sec2day = 1.15741e-5

def calc_bsfc(datadict, W=85., LD=25., Pe=100., rhofuel=750.):
    """ Calculates predicted BSFC

    This method assumes that the datadict is pretrimmed to a stable
    portion of flight.

    INPUTS
    ------
    datadict - dict with flight data from raspberri pi
    W (optional) - aircraft weight during flight segment in pounds
    LD (optional) - lift to drag ratio during flight segement
    Pe (optional) - electrical power draw during flight segment in Watts
    rhofuel (optional) - fuel density in kg/m^3

    OUTPUTS
    -------
    BSFC - predicted bsfc during flight

    """

    W *= lbs2N
    V = np.average(datadict["speed"]["values"])*kts2ms
    totalfuel = datadict["totalfuel"]["values"]
    fuel = (totalfuel[-1] - totalfuel[0])*ml2m3
    fuelrate = datadict["fuelflow"]["values"]*min2sec
    time = datadict["timeelapsed"]["values"]
    totaltime = time[-1] - time[0]
    fuelt = np.trapz(fuelrate, time, np.diff(time)[0])*ml2m3
    assert abs(fuel - fuelt) < 0.001, ("Fuel volume doesn't match; " +
                                       "could be bad data")

    bsfc = fuel*rhofuel/(W/LD*V)/totaltime
    bsfc *= sm2kgkWhr
    return bsfc

def endurance(bsfc, Wtotal=150., Wdry=55., LD=25., CL=0.6, Pe=100., rhofuel=750., S=22.5, rho=0.770816, dt=10):
    """ Calculates predicted endurance

    Time stepping method to calculate endurance

    INPUTS
    ------
    bsfc - BSFC of the aircraft in s^2/m^2
    W (optional) - total aircraft weight in pounds
    LD (optional) - lift to drag ratio assumed constant
    CL (optional) - lift coefficient assumed constant
    Pe (optional) - electrical power draw in Watts
    rhofuel (optional) - fuel density in kg/m^3
    S (optional) - wing area in ft^2
    rho (optional) - air density at flight altitude in kg/m^3
    dt (optional) - time step in seconds

    OUTPUTS
    -------
    BSFC - predicted bsfc during flight

    """

    S *= ft2m**2
    bsfc *= 1./sm2kgkWhr
    Wtotal *= lbs2N
    Wdry *= lbs2N

    W, V, Pa, Ptot, Wfuel, t = [], [], [], [], [], []

    V.append((2*Wtotal/rho/S/CL)**0.5)
    t.append(0)
    Pa.append(Wtotal/LD*V[-1])
    Ptot.append(Pa[-1] + Pe)
    Wfuel.append(Ptot[-1]*bsfc*dt)
    W.append(Wtotal- Wfuel[-1])

    i = 1
    while W[-1] >= Wdry:
        V.append((2*W[-1]/rho/S/CL)**0.5)
        t.append(dt + t[-1])
        Pa.append(W[-1]/LD*V[-1])
        Ptot.append(Pa[-1] + Pe)
        Wfuel.append(Ptot[-1]*bsfc*dt)
        W.append(W[-1]- Wfuel[-1])
        i += 1

    print "Total Endurance: %.1f [days]" % (t[-1]*sec2day)
    return t, W, V, Pa, Ptot

def plot_endurance(time, weight, speed, power_air, power_tot):

    fig, ax = plt.subplots(4)
    ax[0].plot(time, weight)
    ax[0].grid()
    ax[0].set_ylabel("Weight [N]")
    ax[1].plot(time, speed)
    ax[1].grid()
    ax[1].set_ylabel("Speed [m/s]")
    ax[2].plot(time, power_air)
    ax[2].grid()
    ax[2].set_ylabel("Power [W]")
    ax[3].plot(time, power_tot)
    ax[3].grid()
    ax[3].set_ylabel("Tot Power [W]")
    ax[3].set_xlabel("Time [sec]")
    return fig, ax

if __name__ == "__main__":
    D = parse_pidata("20180507191728flightlog-star/2018-05-07-19:17:28-jho_command.log")

    trim_data(D, [3100, 5200])
    print_flightstats(D)
    f, a = plot_params(D, ["speed", "altitude", "pitch", "fuelflow"])
    f.savefig("flight.pdf", bbox_inches="tight")

    f, a = plt.subplots()
    plt.scatter(D["rpm"]["values"], D["fuelflow"]["values"])
    f.savefig("fuelflowrmp.pdf", bbox_inches="tight")

    trim_data(D, [4400, 4600])

    BSFC = calc_bsfc(D)
    T, w, v, pa, ptot = endurance(4.)

    f, a = plot_endurance(T, w, v, pa, ptot)
    f.savefig("endurance.pdf", bbox_incehs="tight")


