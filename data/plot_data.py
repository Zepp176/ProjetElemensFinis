import numpy as np
import matplotlib.pyplot as plt

filename = "couples.txt"
X = np.loadtxt(filename)

t = X[:,0]
theta = X[:,1]
torque = X[:,2]
omega = X[:,3]
zero = np.zeros(len(t))
switch = 3.1415/6
Xone = [switch, switch]
Y = [-0.002, 0.005]
X = [-0.204, -0.204]

fig, ax = plt.subplots(2,1)
ax[0].plot(theta, torque)
ax[0].plot(theta, zero, '--k')
for i in range(6):
    X[0] += Xone[0]
    X[1] += Xone[1]
    ax[0].plot(X, Y, '--k')
ax[0].set_title("Couple")

ax[1].plot(t, omega)

plt.show()
