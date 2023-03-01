import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import stats

data = pd.read_csv('Compression_x_data.csv')
data["y"] = -data["y"]+data["y"][0]
data = data.drop(data[data.y < 0.05].index)

nx = 5
ny = 4
nz = 4

rad = 1.0

dx = 1.75
dy = 1.75
dz = 1.75

Lx = 2*rad + (nx-1)*dx
Ly = 2*rad + (ny-1)*dx
Lz = 2*rad + (nz-1)*dx

A = Lx*Lz 

x = data["y"].to_numpy()
y = data["F"].to_numpy()*Ly/A

res = stats.linregress(x, y)

r2 =res.rvalue**2

plt.title('Compression test on x axis')
plt.plot(x, y, label='Force on wall')
plt.plot(x, res.intercept + res.slope*x, 'r', label=r'Fit: intercept=%5.6g, slope=%5.6g, $r^2$=%5.6g' % tuple([res.intercept,res.slope,r2]))
plt.xlabel('Displacement')
plt.ylabel(r'$F \cdot L/A$')
plt.legend()
plt.show()

