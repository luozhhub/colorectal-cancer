#!/usr/bin/pyhton3
import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def drawSubplot(sample_name, data, location=None):
    fontsize = 12
    ax = plt.subplot(location)
    ax.set_title(sample_name, fontsize=fontsize)
    #ax.set_xlabel(vectorDict["1"][0], fontsize=fontsize)
    #ax.set_ylabel(vectorDict["2"][0], fontsize=fontsize)
    plt.hist(data)


df = pd.read_table("count_table.txt", sep="\t", header="infer")
name_list = df.columns




data = df["823_control_1.txt"][:-5]
print (data[0:5])


plt.figure(1, figsize=(20, 5))
i = 0
for name in name_list[10:15]:
    i += 1
    #i = (i + 1) % 5
    #if i == 0:
    #    i = 5
    data = np.log(df[name][0:-5] + 1)
    drawSubplot(sample_name=name, data=data, location="15%s"%i)
#plt.hist(data, bins=30, density=True)
#plt.hist(data)
#plt.xlim([0, 25000])
#plt.ylim([0, 100000])
#plt.hist(np.log(data + 1))

plt.savefig("result2.png")


