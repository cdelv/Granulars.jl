import pandas as pd 
import matplotlib.pyplot as plt

c = pd.read_csv('data.csv')
julia = pd.read_csv('spinning_top.csv')

plt.plot(c['t'], c['a1'],label='C++')
plt.plot(julia['q3(t)'], julia['q1(a1)'],label='Julia')
plt.legend()
plt.xlabel('time')
plt.ylabel(r'$\theta$')
plt.show()