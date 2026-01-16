# this is for 2A + B = C reaction
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

k1 = 1
k2 = 1


conc0 = np.array([1, 1, 0])

def dfdt(t, conc):
    A, B, C = conc
    r1 = k1*A*A*B
    r2 = k2*C
    dAdt = 2*(-r1 + r2)
    dBdt = -r1 + r2
    dCdt = r1 - r2
    return [dAdt, dBdt, dCdt]

conc = solve_ivp(fun = dfdt, 
                 t_span=(0, 10), 
                 y0 = conc0,
                 )


plt.plot(conc.t, conc.y[0], label='A(t)', color='blue')
plt.plot(conc.t, conc.y[1], label='B(t)', color='red')
plt.plot(conc.t, conc.y[2], label='C(t)', color='green')

plt.xlabel('Time (t)')
plt.ylabel('Values')
plt.title('A+B = C reaction somulation')
plt.legend() 
plt.grid(True)
plt.show()
    