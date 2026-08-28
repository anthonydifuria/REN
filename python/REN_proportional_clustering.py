import numpy as np
import matplotlib.pyplot as plt

# REN parameters
k = 8
lam = 1.6
q = 2.0

# Cluster parameters
theta0 = np.pi / 2   # cluster center
mu = 0.2             # intensity
kappa = 8.0          # width

# Angles
theta = np.linspace(0, 2*np.pi, 4000)
dtheta = theta[1] - theta[0]

# Angular density
bump = np.exp(kappa * np.cos(theta - theta0))
bump /= bump.mean()

w = (1.0 - mu) + mu * bump

# CDF
F = np.cumsum(w) * dtheta

# Normalization
phi = 2*np.pi * (F - F.min()) / (F.max() - F.min())

# Final REN
r = np.exp(-lam * np.abs(np.sin(k * phi))**q)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.plot(theta, r)
ax.set_title("REN proportional cluster")
ax.set_rlim(0, 1)
ax.grid(True)
plt.show()
