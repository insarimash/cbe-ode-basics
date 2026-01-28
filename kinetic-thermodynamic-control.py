# this is for A = B(thermodynamic favorite) A = C(kinetic favorite) reactions
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.linalg import solve

R = 8.314
T = 300
A_preexp = 1000
T_onset = -1.0
A0 = 1
B0 = 0
C0 = 0
T_90 = np.zeros(900)

Ea1b = 60000
Ea2b = 100000

Ea1c = 30000
Ea2c = 50000

# k1b = A*np.exp(-Ea1b/(R*T))
# k2b = A*np.exp(-Ea2b/(R*T))

# k1c = A*np.exp(-Ea1c/(R*T))
# k2c = A*np.exp(-Ea2c/(R*T))


conc0 = np.array([A0, B0, C0])

def dfdt(t, conc, k1b, k2b, k1c, k2c):
    A, B, C = conc
    r1b = k1b*A
    r2b = k2b*B
    r1c = k1c*A
    r2c = k2c*C
    dAdt = - r1b - r1c + r2b + r2c
    dBdt = r1b - r2b
    dCdt = r1c - r2c
    return [dAdt, dBdt, dCdt]

def get_equilibrium(k1b, k2b, k1c, k2c):
    A_array = np.array([[1, 1, 1],
                     [k1b, -k2b, 0],
                     [k1c, 0, -k2c]])
    total = A0 + B0 + C0
    b_array = np.array([total, 0, 0])

    result = solve(A_array, b_array)
    return result
    

# conc = solve_ivp(fun = dfdt, 
#                  t_span=(0, 100), 
#                  y0 = conc0,
#                  )


# plt.plot(conc.t, conc.y[0], label='A(t)', color='blue')
# plt.plot(conc.t, conc.y[1], label='B(t)', color='red')
# plt.plot(conc.t, conc.y[2], label='C(t)', color='green')

# plt.xlabel('Time (t)')
# plt.ylabel('Values')
# plt.title('kinetic-thermodynamic control simulation')
# plt.legend() # This shows the labels we defined in plt.plot
# plt.grid(True)
# plt.show()






conc = np.zeros((900, 3))
T0 = 100


for i in range(900):
    T = T0 + i
    k1b = A_preexp*np.exp(-Ea1b/(R*T))
    k2b = A_preexp*np.exp(-Ea2b/(R*T))

    k1c = A_preexp*np.exp(-Ea1c/(R*T))
    k2c = A_preexp*np.exp(-Ea2c/(R*T))

    A_eq, B_eq, C_eq = get_equilibrium(k1b, k2b, k1c, k2c)

    solution = solve_ivp(fun = dfdt, 
                 t_span=(0, 100), 
                 y0 = conc0,
                 args=(k1b, k2b, k1c, k2c)
                 )
    t = solution.t
    B_sol = solution.y[1]
    C_sol = solution.y[2]
    for j in range(len(t)):
        # if(B_sol[j] >= 0.9 * B_eq   or C_sol[j] >= 0.9 * C_eq):
        #     T_90[i] = t[j]
        #     break
        err = np.linalg.norm(solution.y[:, j] - np.array([A_eq, B_eq, C_eq]))
        err0 = np.linalg.norm(conc0 - np.array([A_eq, B_eq, C_eq]))
        if err <= 0.1 * err0:
            T_90[i] = t[j]
            break
    
    conc[i] = [ solution.y[0][-1], solution.y[1][-1], solution.y[2][-1] ]
    if(T_onset == -1 and float(solution.y[0][-1])/A0 <= 0.99):
        T_onset = T

# print("T_onset = " + str(T_onset))
T_array = np.arange(100, 1000)

# plt.plot(T_array, conc[:, 0], label='A(t)', color='blue')
# plt.plot(T_array, conc[:, 1], label='B(t)', color='red')
# plt.plot(T_array, conc[:, 2], label='C(t)', color='green')

# plt.xlabel('Temperature')
# plt.ylabel('Values')
# plt.title('Final concentrations/Temperature dependancy')
# plt.legend() # This shows the labels we defined in plt.plot
# plt.grid(True)
# plt.show()

plt.plot(T_array, T_90)
plt.xlabel('Temperature')
plt.ylabel('Time to 90')
plt.title("Time to 90%% equilibrium vs. Temperature")
plt.legend()
plt.grid(True)
plt.show()