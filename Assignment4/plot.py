import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000])
y_b = np.array([0, 0.153, 0.486, 0.912, 1.297, 1.768, 2.194, 2.716, 3.187, 3.730, 4.216])
y_old = np.array([0, 0.986, 3.844, 8.778, 15.552, 24.223, 34.669, 47.342, 61.740, 79.292, 96.145])

line1 = plt.plot(x,y_b, label="Barnes hut")
line2 = plt.plot(x, y_old, label="Normal")


plt.xlabel("Number of particles[N]")
plt.ylabel("Time[s]")
plt.title("Complexity graph for both methods")
plt.legend()

plt.show()
