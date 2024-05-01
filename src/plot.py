import sys
import numpy as np
import matplotlib.pyplot as plt


# list column names
lname = ["time", "e-", "H", "OH", "H2O2", "O2", "O2-", "HO2",
         "H2", "H2O", "OH-", "HO2-", "H+", "O-", "O3-", "total"]


# get filename from arg
if len(sys.argv) < 2:
    print("Usage: python plot.py <filename>")
    sys.exit(1)
fname = sys.argv[1]

# read data from file, separated by space, ignoring comments with #
data = np.loadtxt(fname, comments="#")

# make data accessible by molecule name using a dict
d = {}
for i in range(data.shape[1] - 1):
    d[lname[i]] = data[:, i]


fig = plt.figure()
ax = fig.add_subplot(111)

# plot e-
ax.plot(d["time"], d["e-"], label="e-")
ax.plot(d["time"], d["H"], label="H", linestyle="--")
ax.plot(d["time"], d["OH"], label="OH", linestyle="--")
ax.plot(d["time"], d["H2O2"], label="H2O2", linestyle=":")
ax.plot(d["time"], d["O2"], label="O2", linestyle="-.")
ax.plot(d["time"], d["O2-"], label="O2-", linestyle="-")
ax.plot(d["time"], d["HO2"], label="HO2", linestyle="--")
ax.plot(d["time"], d["H2"], label="H2", linestyle="-.")
ax.plot(d["time"], d["H2O"], label="H2O", linestyle=":")
ax.plot(d["time"], d["OH-"], label="OH-", linestyle="--")
ax.plot(d["time"], d["HO2-"], label="HO2-", linestyle="--")
ax.plot(d["time"], d["H+"], label="H+")
ax.plot(d["time"], d["O-"], label="O-", linestyle="-.", color="grey")
ax.plot(d["time"], d["O3-"], label="O3-", linestyle=":")


# log y scale
ax.set_yscale("log")

# limit to y 1e-9 to 1e-3 mmol/l
ax.set_ylim(1e-14, 1e2)
ax.set_ylabel("Concentration [mol/l]")

# set x label
ax.set_xlabel("Time [s]")
# legend position to the right side
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.grid(True)
plt.show()
