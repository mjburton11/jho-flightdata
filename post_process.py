" process flight data from raspberri pi"
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy import stats
import matplotlib.pyplot as plt
from os import listdir, mkdir, rename
from processdata import parse_pidata, plot_params, trim_data, print_flightstats, save_kml

kts2ms = 0.5144444
lbs2N = 4.44822
ml2m3 = 1e-6
min2sec = 1/60.
ft2m = 0.3048
sm2kgkWhr = 1000.*60.*60.
sec2day = 1.15741e-5
g = 9.8 # m/s^2

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
    BSFC - predicted bsfc during flight in kg/kW/hr

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

def endurance(bsfc, Wtotal=150., Wdry=60., LD=25., CL=0.6, Pe=100., rhofuel=750., S=22.5, rho=0.770816, dt=100):
    """ Calculates predicted endurance

    Time stepping method to calculate endurance

    INPUTS
    ------
    bsfc - BSFC of the aircraft in kg/kW/hr
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
    dict {time          - flight time in days
          weight        - aircraft weight in lbs
          speed         - speed in knots
          flight power  - power required to fly in kW
          total power   - flight power plus avionics power in kW}

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
        Wfuel.append(Ptot[-1]*bsfc*dt*g)
        W.append(W[-1]- Wfuel[-1])
        i += 1

    t, W, V, Pa = np.array(t), np.array(W), np.array(V), np.array(Pa)
    Ptot = np.array(Ptot)

    print "Total Endurance: %.1f [days]" % (t[-1]*sec2day)
    return {"time": {"units": "days", "values": t/60./60./24.},
            "weight": {"units": "lbs", "values": W/lbs2N},
            "speed": {"units": "kts", "values": V/kts2ms},
            "flight power": {"units": "kW", "values": Pa/1000.},
            "total power": {"units": "kW", "values": Ptot/1000.}}

def calc_ld(data, W=85., plotfits=True):

    W *= lbs2N
    m = W/g
    t = data["timeelapsed"]["values"] - data["timeelapsed"]["values"][0]
    ti = np.linspace(t[0], t[-1], 100)
    dt = np.diff(ti)

    h = data["altitude"]["values"]*ft2m
    m, b, rval, _, _ = stats.linregress(t, h)
    hi = m*ti + b
    dh = hi[:-1] - hi[1:]
    if plotfits:
        fig, ax = plt.subplots()
        ax.plot(ti, hi)
        ax.grid()
        ax.scatter(t, h, facecolor="None")
        ax.set_xlabel("Time [sec]")
        ax.set_ylabel("Altitude [m]")
        ax.set_title("R-squared value: %.3f" % (rval**2))
        fig.savefig("htest.pdf")

    r = data["roll"]["values"]
    fr = interp1d(t, r, "cubic")
    ri = fr(ti)
    from savitzky_golay import savitzky_golay
    ri = savitzky_golay(ri, window_size=41, order=2)
    if plotfits:
        fig, ax = plt.subplots()
        ax.plot(ti, ri)
        ax.grid()
        ax.scatter(t, r, facecolor="None")
        ax.set_xlabel("Time [sec]")
        ax.set_ylabel("Roll [radians]")
        fig.savefig("rtest.pdf")
    ri = 0.5*(ri[1:] + ri[:-1])
    cosr = np.cos(ri)

    V = np.average(data["speed"]["values"])*kts2ms

    LD = V/cosr/dh*dt
    if plotfits:
        fig, ax = plt.subplots()
        ax.plot(ti[:-1], LD)
        ax.grid()
        ax.set_xlabel("Time [sec]")
        ax.set_ylabel("L/D")
        fig.savefig("LD_constantV.pdf")

    V = data["speed"]["values"]*kts2ms
    fv = interp1d(t, V, "cubic")
    Vi = fv(ti)
    Vi = savitzky_golay(Vi, window_size=31, order=4)
    if plotfits:
        fig, ax = plt.subplots()
        ax.plot(ti, Vi)
        ax.grid()
        ax.scatter(t, V, facecolor="None")
        ax.set_xlabel("Time [sec]")
        ax.set_ylabel("Speed [m/s]")
        fig.savefig("Vtest.pdf")

    dKE = 0.5*m*(Vi[1:]**2 - Vi[:-1]**2)
    V = 0.5*(Vi[1:] + Vi[:-1])
    LD = 1./((m*g*dh/dt + dKE)/m/g*cosr/V)
    if plotfits:
        fig, ax = plt.subplots()
        ax.plot(ti[:-1], LD)
        ax.grid()
        ax.set_xlabel("Time [sec]")
        ax.set_ylabel("L/D")
        # ax.set_ylim([0, 70])
        fig.savefig("LD.pdf")
    return LD

def plot_endurance(datadict):

    time = datadict["time"]["values"]
    fig, ax = plt.subplots(4)
    ax[3].set_xlabel("Time [" + datadict["time"]["units"] + "]")
    datadict.pop("time")
    for d, x in zip(datadict, ax):
        x.plot(time, datadict[d]["values"])
        x.grid()
        x.set_ylabel(d + " [" + datadict[d]["units"] + "]")
    return fig, ax

if __name__ == "__main__":
    D = parse_pidata("test_12_flightlog/test_12_jho_command.log")

    # trim_data(D, [3367, 3380], tempmax=145)
    # trim_data(D, [3367, 3390], tempmax=145)
    # ld = calc_ld(D)
    # print_flightstats(D)
    # f, a = plot_params(D, ["speed", "altitude", "roll", "heading"])
    # f.savefig("test_12_flightlog/flight.pdf", bbox_inches="tight")

    # f, a = plot_params(D, ["rpm", ["cht1", "cht2"], "pressure"])
    # f.savefig("test_12_flightlog/engine.pdf", bbox_inches="tight")

    # trim_data(D, [1240, 1265], tempmax=145)
    trim_data(D, [600, 3500], tempmax=145)
    # print_flightstats(D)
    # save_kml(D, "test_12_flightlog/test_12.kml")
    f, a = plot_params(D, ["speed", "altitude", "roll", "heading"])
    f.savefig("test_12_flightlog/flight.pdf", bbox_inches="tight")

    f, a = plot_params(D, ["rpm", ["cht1", "cht2"], "pressure"])
    f.savefig("test_12_flightlog/engine.pdf", bbox_inches="tight")

    # trim_data(D, [1500, 1750], tempmax=145)
    # BSFC = calc_bsfc(D, W=82.)
    # enddict = endurance(BSFC)
    # f, a = plot_endurance(enddict)
    # f.savefig("test_12_flightlog/endurance.pdf")

    # trim_data(D, [3300, 3450], tempmax=145)
    # ld = calc_ld(D)
