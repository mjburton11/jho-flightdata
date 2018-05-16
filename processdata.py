" process flight data from raspberri pi"
import pandas as pd
import numpy as np
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
        if "*" in d or "#" in d:
            continue
        if "start" in d or "Connect" in d:
            continue
        st = d.replace("[", "").replace("]", "").split(",")
        nums = [float(s) for s in st]
        if any(nums) == 0.0:
            continue
        dt = datetime.strptime(t, "%Y-%m-%d %H:%M:%S")
        tmel = (dt - data[0][-1]).total_seconds() if data else 0.0
        data.append(nums + [tmel, dt])

    data = np.array(data)

    ms2kts = 1.94384
    ddict = {"altitude": {"units": "ft", "values": -data[:, 0]},
             "heading": {"units": "radians", "values": data[:, 1]},
             "speed": {"units": "kts", "values": data[:, 2]*ms2kts},
             "pitch": {"units": "radians", "values": data[:, 3]},
             "roll": {"units": "radians", "values": data[:, 4]},
             "pressure": {"units": "bar", "values": data[:, 5]},
             "rpm": {"units": "RPM", "values": data[:, 6]},
             "ecuvoltage": {"units": "V", "values": data[:, 7]},
             "cht1": {"units": "C", "values": data[:, 8]},
             "cht2": {"units": "C", "values": data[:, 9]},
             "fuelflow": {"units": "ml/min", "values": data[:, 10]},
             "totalfuel": {"units": "ml", "values": data[:, 11]},
             "gpse": {"units": "degrees", "values": data[:, 12]},
             "gpsn": {"units": "degrees", "values": data[:, 13]},
             "voltage": {"units": "V", "values": data[:, 14]},
             "servovolt": {"units": "V", "values": data[:, 15]},
             "fuelpresence": {"units": "-", "values": data[:, 16]},
             "payloadtemp": {"units": "C", "values": data[:, 17]},
             "mptemp": {"units": "C", "values": data[:, 18]},
             "timeelapsed": {"units": "sec", "values": data[:, 20]},
             "time": {"units": "-", "values": data[:, 21]}
            }

    return ddict

def save_kml(datadict, kmlfile):
    " print out kml file of lat long "
    dfd = {"longitude": datadict["gpse"]["values"],
           "latitude": datadict["gpsn"]["values"]}
    df = pd.DataFrame(dfd)
    kml = simplekml.Kml()
    df.apply(lambda X: kml.newpoint(
        name="1", coords=[(X["longitude"], X["latitude"])]), axis=1)
    kml.save(path=kmlfile)

def print_flightstats(datadict):
    " print import flight characteristics "

    print "Log time: %d [min]" % (datadict["timeelapsed"]["values"][-1]/60.)
    print "Max Altitude: %d [%s]" % (max(datadict["altitude"]["values"]),
                                     datadict["altitude"]["units"])
    print "Max Speed: %.2f [%s]" % (max(datadict["speed"]["values"]),
                                    datadict["speed"]["units"])
    print "Max RPM: %d" % max(datadict["rpm"]["values"])
    chts = [max(datadict["cht1"]["values"]), max(datadict["cht2"]["values"])]
    print "Max CHT: %.2f [%s]" % (max(chts), datadict["cht1"]["units"])

def check_units(datadict, params):
    " check if units are the same "
    assert len(params) == 2, "Just give me 2 params!"
    assert datadict[params[0]]["units"] == datadict[params[1]]["units"], (
        "Units must be the same for double plotting!")

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
            ax[i].legend(params[i])
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
        save_kml(D, dirname + "/coords.kml")
        print_flightstats(D)
        f, a = plot_params(D, ["rpm", ["cht1", "cht2"], "ecuvoltage",
                               "pressure"])
        f.savefig(dirname + "/engine.pdf", bbox_inches="tight")
        f, a = plot_params(D, ["speed", "altitude", "pitch", "fuelflow"])
        f.savefig(dirname + "/flight.pdf", bbox_inches="tight")
        f, a = plot_params(D, ["voltage", "payloadtemp", "mptemp"])
        f.savefig(dirname + "/avionics.pdf", bbox_inches="tight")

