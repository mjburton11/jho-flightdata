" process flight data from raspberri pi"
import pandas as pd
import numpy as np
from numpy import array as ay
from numpy import float32 as flt
import simplekml
from datetime import datetime
import matplotlib.pyplot as plt
from os import listdir, mkdir, rename

def parse_pidata(pifile):
    " parse raspberrypi data file from JHO"
    df = pd.read_csv(pifile, delimiter=";")

    for i, d in enumerate(df[df.columns[2]]):
        if "[" in d:
            istart = i
            break

    data = []
    for d, t in zip(df.iloc[istart:, 2], df.iloc[istart:, 0]):
        if isinstance(d, float) and np.isnan(d):
            continue
        if "*" in d or "#" in d or "L" in d:
            continue
        if "start" in d or "Connect" in d:
            continue
        st = d.replace("[", "").replace("]", "").split(",")
        nums = []
        for s in st:
            if "None" not in s:
                nums.append(float(s))
            else:
                nums.append(np.nan)
        if any(nums) == 0.0:
            continue
        dt = datetime.strptime(t, "%Y-%m-%d %H:%M:%S")
        tmel = (dt - data[0][-1]).total_seconds() if data else 0.0
        data.append(nums + [tmel, dt])

    data = np.array(data)

    fs2kts = 0.592484
    ft2m = 0.3048
    ddict = {"altitude": {"units": "ft", "values": ay(-data[:, 0], dtype=flt)},
             "heading": {"units": "radians",
                         "values": ay(data[:, 1], dtype=flt)},
             "speed": {"units": "kts",
                       "values": ay(data[:, 2]*fs2kts, dtype=flt)},
             "pitch": {"units": "radians", "values": ay(data[:, 3], dtype=flt)},
             "roll": {"units": "radians", "values": ay(data[:, 4], dtype=flt)},
             "pressure": {"units": "bar", "values": ay(data[:, 5], dtype=flt)},
             "rpm": {"units": "RPM", "values": ay(data[:, 6], dtype=flt)},
             "ecuvoltage": {"units": "V", "values": ay(data[:, 7], dtype=flt)},
             "cht1": {"units": "C", "values": ay(data[:, 8], dtype=flt)},
             "cht2": {"units": "C", "values": ay(data[:, 9], dtype=flt)},
             "fuelflow": {"units": "ml/min",
                          "values": ay(data[:, 10], dtype=flt)},
             "totalfuel": {"units": "ml", "values": ay(data[:, 11], dtype=flt)},
             "gpse": {"units": "degrees", "values": ay(data[:, 12], dtype=flt)},
             "gpsn": {"units": "degrees", "values": ay(data[:, 13], dtype=flt)},
             "voltage": {"units": "V", "values": ay(data[:, 14], dtype=flt)},
             "servovolt": {"units": "V", "values": ay(data[:, 15], dtype=flt)},
             "fuelpresence": {"units": "-",
                              "values": ay(data[:, 16], dtype=flt)},
             "payloadtemp": {"units": "C",
                             "values": ay(data[:, 17], dtype=flt)},
             "mptemp": {"units": "C", "values": ay(data[:, 18], dtype=flt)},
             "zaccel": {"units": "m/s^2",
                        "values": ay(data[:, 19]*ft2m, dtype=flt)},
             "xaccel": {"units": "m/s^2",
                        "values": ay(data[:, 20]*ft2m, dtype=flt)},
             "yaccel": {"units": "m/s^2",
                        "values": ay(data[:, 21]*ft2m, dtype=flt)},
             "timeelapsed": {"units": "sec",
                             "values": ay(data[:, 22], dtype=flt)},
             "time": {"units": "-", "values": data[:, 23]}
            }

    return ddict

def save_kml(datadict, kmlfile):
    " print out kml file of lat long "
    dfd = {"longitude": datadict["gpse"]["values"],
           "latitude": datadict["gpsn"]["values"]}
    df = pd.DataFrame(dfd)
    df = df[df.longitude > -65]
    df = df[df.longitude < -75]
    df = df[df.latitude < 35]
    df = df[df.latitude > 45]
    kml = simplekml.Kml()
    df.apply(lambda X: kml.newpoint(
        name="1", coords=[(X["longitude"], X["latitude"])]), axis=1)
    kml.save(path=kmlfile)

def print_flightstats(datadict):
    " print import flight characteristics "

    print "Log time: %d [min]" % ((datadict["timeelapsed"]["values"][-1]
                                   - datadict["timeelapsed"]["values"][0])/60.)
    print "Max Altitude: %d [%s]" % (max(datadict["altitude"]["values"]),
                                     datadict["altitude"]["units"])
    print "Max Speed: %.2f [%s]" % (max(datadict["speed"]["values"]),
                                    datadict["speed"]["units"])
    print "Max RPM: %d" % max(datadict["rpm"]["values"])
    chts = [max(datadict["cht1"]["values"]), max(datadict["cht2"]["values"])]
    print "Max CHT: %.2f [%s]" % (max(chts), datadict["cht1"]["units"])
    ml2gl = 0.000264172
    totalfuel = datadict["totalfuel"]["values"]
    fuel = (totalfuel[-1] - totalfuel[0])*ml2gl
    print "Total fuel burn: %.2f [gallon]" % fuel

def check_units(datadict, params):
    " check if units are the same "
    assert len(params) == 2, "Just give me 2 params!"
    assert datadict[params[0]]["units"] == datadict[params[1]]["units"], (
        "Units must be the same for double plotting!")

def trim_data(datadict, trange, rpmmax=9000, tempmax=170, barmax=5,
              pitchmin=-50, pitchmax=50):
    "trim flight data to specific range"

    assert len(trange) == 2, "Specify upper and lower range values"
    for i in range(2):
        t = datadict["timeelapsed"]["values"]
        ind = t > trange[i] if i == 0 else t < trange[i]
        for d in datadict:
            if d == "timeelapsed":
                continue
            datadict[d]["values"] = datadict[d]["values"][ind]
        datadict["timeelapsed"]["values"] = t[ind]

    datadict["rpm"]["values"][datadict["rpm"]["values"] > rpmmax] = np.nan
    datadict["cht1"]["values"][datadict["cht1"]["values"] > tempmax] = np.nan
    datadict["cht2"]["values"][datadict["cht2"]["values"] > tempmax] = np.nan
    datadict["pressure"]["values"][
        datadict["pressure"]["values"] > barmax] = np.nan
    datadict["pitch"]["values"][datadict["pitch"]["values"] < pitchmin] = np.nan
    datadict["pitch"]["values"][datadict["pitch"]["values"] > pitchmax] = np.nan

def plot_params(datadict, params):
    " plot flight parameters and return figure "

    N = len(params)
    fig, ax = plt.subplots(N)
    for i in range(N):
        x = datadict["timeelapsed"]["values"]
        if isinstance(params[i], list):
            check_units(datadict, [params[i][0], params[i][1]])
            yunits = datadict[params[i][0]]["units"]
            for j in range(len(params[i])):
                y = datadict[params[i][j]]["values"]
                ax[i].plot(x, y)
            ax[i].set_ylabel("[%s]" % yunits)
            # ax[i].legend(params[i])
        else:
            y = datadict[params[i]]["values"]
            yunits = datadict[params[i]]["units"]
            ax[i].plot(x, y)
            ax[i].set_ylabel("%s [%s]" % (params[i], yunits))
        ax[i].grid()
    ax[-1].set_xlabel("time [%s]" % datadict["timeelapsed"]["units"])
    return fig, ax

if __name__ == "__main__":
    files = listdir("./")
    for f in files:
        if not "jho_command" in f:
            continue
        dirname = f.replace("jho_command.log", "flightlog").replace(
            "-", "").replace(":", "")
        mkdir("./" + dirname)
        rename(f, dirname + "/" + f)
        D = parse_pidata(dirname + "/" + f)
        f, a = plot_params(D, ["rpm", ["cht1", "cht2"], "ecuvoltage",
                               "pressure"])
        f.savefig(dirname + "/engine.pdf", bbox_inches="tight")
        f, a = plot_params(D, ["speed", "altitude", "pitch", "fuelflow"])
        f.savefig(dirname + "/flight.pdf", bbox_inches="tight")
        f, a = plot_params(D, ["voltage", "payloadtemp", "mptemp"])
        f.savefig(dirname + "/avionics.pdf", bbox_inches="tight")

