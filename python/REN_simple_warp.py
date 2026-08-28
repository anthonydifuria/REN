import numpy as np
import matplotlib.pyplot as plt

# REN parameters
k = 8
lam = 1.6
q = 2.0

# Deformation parameters
theta0 = np.pi / 2   # cluster center (90 degrees)
mu = 0.6             # deformation intensity

# Angles
theta = np.linspace(0, 2*np.pi, 4000)

# Simple warp
phi = theta + mu * np.sin(theta - theta0)

# Warped REN
r = np.exp(-lam * np.abs(np.sin(k * phi))**q)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.plot(theta, r)
ax.set_title("REN simple warp")
ax.set_rlim(0, 1)
ax.grid(True)
plt.show()
