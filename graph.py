import matplotlib.pyplot as plt
import numpy as np
import csv
import os

directory = '../src/summaries'
summaries = os.listdir(directory)

for summary in summaries:
    averages = []
    stddev = []

    if "averages.csv" in summary:
        with open("/".join([directory, summary]), newline='') as file:
            reader = csv.reader(file, delimiter=' ')
            for row in reader:
                averages.append(float(row[0]))
                stddev.append(float(row[1]))

        fig, ax = plt.subplots()

        ax.plot(np.arange(len(averages)) + 1, averages)
        ax.errorbar(np.arange(len(averages)) + 1, averages, stddev, linestyle='None', marker='None')
        ax.set_title("Average Average Fitness:" + summary)
        ax.set_ylabel("Fitness")
        ax.set_xlabel("Generation")

        plt.savefig(summary.split('.csv')[0] + '.png')
