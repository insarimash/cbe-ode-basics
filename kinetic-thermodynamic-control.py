# this is for A = B(thermodynamic favorite) A = C(kinetic favorite) reactions
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
R = 8.314
T = 300
A = 1000

Ea1b = 60000
Ea2b = 100000

Ea1c = 30000
Ea2c = 50000

k1b = A*np.exp(-Ea1b/(R*T))
k2b = A*np.exp(-Ea2b/(R*T))

k1c = A*np.exp(-Ea1c/(R*T))
k2c = A*np.exp(-Ea2c/(R*T))


conc0 = np.array([1, 0, 0])

def dfdt(t, conc):
    A, B, C = conc
    r1b = k1b*A
    r2b = k2b*B
    r1c = k1c*A
    r2c = k2c*C
    dAdt = - r1b - r1c + r2b + r2c
    dBdt = r1b - r2b
    dCdt = r1c - r2c
    return [dAdt, dBdt, dCdt]

conc = solve_ivp(fun = dfdt, 
                 t_span=(0, 100), 
                 y0 = conc0,
                 )


plt.plot(conc.t, conc.y[0], label='A(t)', color='blue')
plt.plot(conc.t, conc.y[1], label='B(t)', color='red')
plt.plot(conc.t, conc.y[2], label='C(t)', color='green')

plt.xlabel('Time (t)')
plt.ylabel('Values')
plt.title('kinetic-thermodynamic control simulation')
plt.legend() # This shows the labels we defined in plt.plot
plt.grid(True)
plt.show()
    