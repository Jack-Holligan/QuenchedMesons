from argparse import ArgumentParser
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import random
from scipy import interpolate
import yaml


plt.style.use("styles/paperdraft.mplstyle")

parser = ArgumentParser()
parser.add_argument("--beta", type=float)
parser.add_argument("--metadata", default=None)
parser.add_argument("--bootstrap1", type=int, default=200)
parser.add_argument("--bootstrap2", type=int, default=50)
parser.add_argument("--flow_file", default=None)
parser.add_argument("--output_file_main", default=None)
parser.add_argument("--output_file_plot", default=None)
args = parser.parse_args()

flowtime_values = []
t2Eplaq_values = []
t2Esym_values = []
with open(args.flow_file, "r") as f:
    for line in f:
        if line.startswith("[SYSTEM][0]Gauge group: SP"):
            Nc = int(line.split()[-1].strip("SP()"))
            scale = 0.28125 * (Nc + 1) / 5

        if line.startswith("[GEOMETRY_INIT][0]Global size is"):
            num_sites_t, num_sites_x, num_sites_y, num_sites_z = map(int, line.split()[-1].split("x"))
            if not (num_sites_x == num_sites_y == num_sites_z):
                raise ValueError("We asssume NX == NY == NZ")

        if line.startswith("[MAIN][0]WF number of measures:"):
            num_flowtimes = int(line.split()[4]) + 1

        if line.startswith("[WILSONFLOW][0]WF (t,E,t2*E,Esym,t2*Esym,TC) ="):
            split_line = line.split()
            if (flowtime := float(split_line[3])) > max(flowtime_values + [-1]) and len(flowtime_values) < num_flowtimes:
                flowtime_values.append(flowtime)
            if split_line[3] == "0.0000000000000000e+00":
                if len(t2Eplaq_values) > 0 and len(t2Eplaq_values[-1]) != len(flowtime_values):
                    breakpoint()
                    raise ValueError("Mismatched flow lengths.")
                t2Eplaq_values.append([])
                t2Esym_values.append([])
            elif len(t2Eplaq_values[-1]) >= num_flowtimes:
                continue
            t2Eplaq_values[-1].append(float(split_line[5]))
            t2Esym_values[-1].append(float(split_line[7]))

flowtimes = np.asarray(flowtime_values)
t2Eplaq = np.asarray(t2Eplaq_values)
t2Esym = np.asarray(t2Esym_values)

cnfg, nmeas = t2Eplaq.shape

with open(args.metadata, "r") as f:
    metadata = yaml.safe_load(f)
random.seed(metadata[f"Sp{Nc}_S{num_sites_x}T{num_sites_t}B{args.beta}"]["seed"])

def tInterval(EsymmAve, scale):
    """Returns the ARRAY NUMBERS of flow time between which the scale lives"""
    i = 0
    while EsymmAve[i] < scale:
        i += 1
    return [i - 1, i]


def bootstrap(boots, t2E):
    """Selects "boots" configurations at random then averages t**2E at equal flow times"""
    bootarray = []
    t2E_samples = np.vstack([t2E[random.randint(0, cnfg - 1), :] for _ in range(boots)])
    return t2E_samples.mean(axis=0)


def w(EsymmAve, time, i):
    """Computes the value of W(t) for a given array number, flowtime array and t2Esymm array"""
    diff = EsymmAve[i] - EsymmAve[i - 1]
    diff /= time[i] - time[i - 1]
    diff *= time[i - 1]
    return diff


def wscale(EsymmAve, time):
    """Makes an array of W values for given t2Esymm and time arrays."""
    values = np.zeros(len(EsymmAve) - 1)
    for i in range(1, len(values) + 1):
        values[i - 1] = w(EsymmAve, time, i)
    return values


def flowtime(EsymmAve, time, scale):
    """Computes the flow time that corresponds to a given scale"""
    window = tInterval(EsymmAve, scale)
    t1 = time[window[0]]
    t2 = time[window[1]]
    E1 = EsymmAve[window[0]]
    E2 = EsymmAve[window[1]]
    x = [t1, t2]
    y = [E1, E2]
    f = interpolate.interp1d(
        y, x
    )  # x and y are reversed since we want to find the x-value for a given y-value.
    return f(scale)


scalevalues = np.zeros(args.bootstrap2)

for n in range(args.bootstrap2):
    EsymmAve = bootstrap(args.bootstrap1, t2Esym)
    wilson = wscale(EsymmAve, flowtimes)
    tvalue = flowtime(wilson, flowtimes, scale)
    scalevalues[n] = tvalue

mu1, std1 = norm.fit(np.sqrt(scalevalues))
print("Clover: %f %f" % (mu1, std1))

if args.output_file_main:
    with open(args.output_file_main, "w") as f:
        f.write("---------\n")
        f.write("beta %f\nw %f\ndw %f\nW0 %f\n" % (args.beta, mu1, std1, scale))
        f.write("---------\n")

average1 = t2Esym.mean(axis=0)
error1 = t2Esym.std(axis=0)

wvalues1 = []
for i in range(1, cnfg + 1):
    wvalue = w(average1, flowtimes, i)
    wvalues1.append(wvalue)

tw1 = flowtimes[:-1]

for n in range(args.bootstrap2):
    EAve = bootstrap(args.bootstrap1, t2Eplaq)
    wilson = wscale(EAve, flowtimes)
    tvalue = flowtime(wilson, flowtimes, scale)
    scalevalues[n] = tvalue

mu2, std2 = norm.fit(np.sqrt(scalevalues))
print("Plaquette: %f %f" % (mu2, std2))

average2 = t2Eplaq.mean(axis=0)
error2 = t2Eplaq.std(axis=0)

wvalues2 = [w(average2, flowtimes, i) for i in range(1, cnfg + 1)]

tw2 = flowtimes[:-1]

plt.errorbar(tw2, wvalues2, marker="o", label="Plaquette", color="green")
plt.errorbar(tw1, wvalues1, marker="^", label="Clover", color="blue")
plt.ylabel(r"$\mathcal{W}(t)$")
plt.xlabel(r"$t$")
plt.xlim(0, 10.0)
plt.ylim(0, 1.2 * wvalues1[10])
plt.vlines(x=mu1**2, color="blue", ymin=0, ymax=scale)
plt.hlines(y=scale, color="blue", xmin=0, xmax=mu1**2)
plt.vlines(x=mu2**2, color="green", ymin=0, ymax=scale)
plt.hlines(y=scale, linestyle="-", color="green", xmin=0, xmax=mu2**2)
plt.title(r"$Sp(%d), %d^3\times%d, \beta=%1.1f$" % (Nc, num_sites_x, num_sites_t, args.beta))
plt.grid()
plt.legend(loc="best")
plt.savefig(args.output_file_plot)
plt.close()
