import matplotlib.pyplot as plt
import numpy as np
import csv

xs = []
ys = []

with open('../src/best_fit.csv', newline='') as file:
    reader = csv.reader(file, delimiter=' ')
    for row in reader:
        xs.append(float(row[0]))
        ys.append(float(row[1]))

fig, ax = plt.subplots()


ax.scatter(xs, ys)

t = np.linspace(0, np.pi * 2, 100)

ax.plot(np.cos(t), np.sin(t), linewidth=1)
plt.gca().set_aspect('equal')

fig.savefig('best_fit.png')
