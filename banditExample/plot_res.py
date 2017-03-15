import cPickle as pickle

import matplotlib.pyplot as plt

import numpy as np


with open("all_res.p", "r") as f:
    all_res = pickle.load(f)


all_res_handle = [np.mean(i) for i, j in all_res]
all_res_value = [np.mean(j) for i, j in all_res]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot([2 * (i + 1) for i in range(50)],
         all_res_handle[:-1], c = "black")
ax.axhline(y = all_res_handle[-1], c = "red")
fig.savefig("correct_pulls.png")

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot([2 * (i + 1) for i in range(50)],
         all_res_value[:-1], c = "black")
ax.axhline(y = all_res_value[-1], c = "red")
fig.savefig("estimated_value.png")
