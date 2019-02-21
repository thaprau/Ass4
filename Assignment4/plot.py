import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200])
y = np.array([0, 0.242, 0.965, 2.12, 3.625, 5.463, 7.360, 9.822, 12.307, 14.786, 17.368, 20.285])

plt.plot(x,y)
plt.show()
