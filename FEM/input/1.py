import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def is_well_point(point, wells):
    for well in wells:
        if ((point['x'] == well['x0'] or point['x'] == well['x1']) and (point['y'] == well['y0'] or point['y'] == well['y1'])):
            return True
    return False

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
bc1_nodes = []
grid = []
wells = []
with open("BC1.txt", "r") as f:
    next(f)
    for line in f:
        for item in line.split():
            if item.isdecimal():
                bc1_nodes.append(int(item))
with open("grid.txt", 'r') as f:
    next(f)
    for line in f:
        if line == "\n": continue
        lf = line.split(" ")
        x = lf[0]
        y = lf[1]
        grid.append({"x" : float(x), "y" : float(y)})
with open("wells.txt", "r") as f:
    for line in f.readlines()[1:3]:
        items = line.split()
        wells.append({"x0" : float(items[0]), "x1" : float(items[1]), "y0" : float(items[2]), "y1" : float(items[3])})

for ind, point in enumerate(grid):
    if ind in bc1_nodes:
        plt.plot(point['x'], point['y'], "rx")
    elif not is_well_point(point, wells):
        plt.plot(point['x'], point['y'], "ko")

for well in wells:
        plt.plot(well["x0"], well["y0"], "gP")
        plt.plot(well["x0"], well["y1"], "gP")
        plt.plot(well["x1"], well["y0"], "gP")
        plt.plot(well["x1"], well["y1"], "gP")

plt.show()        