import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Import Data
# ----------------------------
T = 341.0 # K 
P = 1700.0  # kPa

z = np.array([0.61, 0.28,0.11])

# Critical properties
Pc = np.array ([4600, 3800, 2108]) 
Tc = np.array ([190.6, 425.2, 617.7]) 
omega = np.array ([0.008, 0.193, 0.49])


# ----------------------------
# Wilson K-values
# ----------------------------
def wilson_K(T, P, Tc, Pc, omega):
    K = []
  
    for i in range(len(Tc)):
      K_i = (Pc[i] /P) * np.exp(5.373 * (1 + omega[i]) * (1 - Tc[i] /T))
      K.append(K_i)
 
  return np.array(K)

# ----------------------------
# Rachord-Rice Function
# ----------------------------
def rr_function(beta, z, K):
    total = 0 

    for i in range(len(z)):
        total += (z[i] * (K[i] - 1)) / (1 + beta * (K[i] - 1))

    return total

# ----------------------------
# Bisection Method
# ----------------------------
def solve_rr(z, K):

    beta_low = 0.0
    beta_high = 1.0

    for _ in range(100):

        beta_mid = (beta_low + beta_high) / 2
        f_mid = rr_function(beta_mid, z, K)

        if abs(f_mid) < 1e-6:
            return beta_mid

        if f_mid > 0:
            beta_high = beta_mid
        else: 
            beta_low = beta_mid
          
    return beta_mid
# ----------------------------
# Main Calculation
# ----------------------------
K = wilson_K(T, P, Tc, Pc, omega)

beta_v = solve_rr(z, K)

# Phase compositions
y = []
x = []

for i in range(len(z)):
    y_i = (K[i] * z[i]) / ( 1 + beta_v * (K[i] - 1))
    x_i = y_i / K[i]

    y.append(y_i)
    x.append(x_i)

y = np.array(y)
x = np.array(x)

# Phase Amounts
V = beta_v
  
            
L = 1 - V

# ----------------------------
# Output
# ----------------------------   
print("Beta V", round(V, 2))
print("Beta L", round(L, 2))

print("Liquid Composition", np.round(x, 2))
print("Vapor Composition", np.round(y, 2))

# ----------------------------
# Plot
# ----------------------------
beta_vals = np.linspace(0, 1, 100)
f_vals = []

for b in beta_vals:
    f_vals.append(rr_function(b, z, K))

plt.plot(beta_vals, f_vals)
plt.axhline(0)
plt.xlabel("Beta V")
plt.ylabel("RR Function")
plt.title( "RR Function v Beta V")
plt.grid()
plt.show()
