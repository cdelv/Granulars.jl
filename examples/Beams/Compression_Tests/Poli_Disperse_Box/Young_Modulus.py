import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import stats

data = pd.read_csv('Compression_z_data.csv')
data["y"] = -data["y"]+data["y"][0]
data = data.drop(data[data.y < 0.1].index)
data = data.drop(data[data.y > 0.116].index)

# z ->
Lx = 6.88
Ly = 6.8
Lz = 6.4

# x ->
#Lx = 6.4
#Ly = 6.8
#Lz = 6.88

A = Lx*Lz 

x = data["y"].to_numpy()
y = data["F"].to_numpy()*Ly/A

res = stats.linregress(x, y)
r2 =res.rvalue**2

plt.title('Compression test on z axis')
plt.plot(x, y, label='Force on wall')
plt.plot(x, res.intercept + res.slope*x, 'r', label=r'Fit: intercept=%5.6g, slope=%5.6g, $r^2$=%5.6g' % tuple([res.intercept,res.slope,r2]))
plt.xlabel('Displacement')
plt.ylabel(r'$F \cdot L/A$')
plt.legend()
plt.show()

