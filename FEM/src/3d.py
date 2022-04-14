import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np

x=[]
y=[]
z=[]

with open("result.txt", "r") as f:
    for line in f.readlines()[3:-2]:
        if (line.startswith("|   ||P-T||/||T||   |") or line.startswith("+--------------------")): continue
        arr = list(line.split("|")[1:])
        arr = [elem.strip() for elem in arr]
        x.append(float(arr[0]))
        y.append(float(arr[1]))
        z.append(float(arr[2]))

x = np.asarray(x)
y = np.asarray(y)
z = np.asarray(z)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x, y, z, linewidth=0, antialiased=False, label='pressure', shade=True)
plt.show()
#, cmap=cm.terrain