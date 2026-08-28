import numpy as np
import matplotlib.pyplot as plt

# REN parameters
k = 8
lam = 1.6
q = 2.0

# Angles
theta = np.linspace(0, 2*np.pi, 4000)

# REN formula
r = np.exp(-lam * np.abs(np.sin(k * theta))**q)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.plot(theta, r)
ax.set_title("REN base")
ax.set_rlim(0, 1)
ax.grid(True)
plt.show()
