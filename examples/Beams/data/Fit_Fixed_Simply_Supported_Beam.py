import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# https://calcresource.com/statics-simple-beam-diagrams.html
# https://mechanicalc.com/reference/beam-analysis

n = 50
rad = 0.5
m = 0.01
g = -10.0

E = 800000.0
I = 0.012767628893730067
L = (2*rad - 0.6*rad)*n

w = m*g*n/L

a = w/(24*E*I)
b = -2*L*w/(24*E*I)
c = w*L**2/(24*E*I)

print("a = ", a)
print("b = ", b)
print("c = ", c)

def func(x, a, b, c):
    return a*x**4 + b*x**3 + c*x**2

data = pd.read_csv('Data_Fixed_Simply_Supported_Beam.csv')
data['X'] -= data['X'][0]
data['Y'] -= data['Y'][0]

X = data['X'].to_numpy()
Y = data['Y'].to_numpy()

popt, pcov = curve_fit(func, X, Y)
print(popt)

plt.title('Fixed Simply Supported Beam')
plt.xlabel('x')
plt.ylabel('Deflection')
plt.plot(X, Y, 'o', color='green', label='Data', alpha=0.7)
plt.plot(X, Y, 'o', color='green', label=r'$\eta(x) = ax^4 + bx^3 + cx^2$', alpha=0.0)
plt.plot(X, func(X, *popt), 'r-', label='fit: a=%5.6g, b=%5.6g, c=%5.6g' % tuple(popt))
plt.plot(X, func(X, a, b, c), 'b-', label='sol: a=%5.6g, b=%5.6g, c=%5.6g' % tuple([a,b,c]))

plt.legend()
plt.show()