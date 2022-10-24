#!/usr/bin/python3

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    nsteps = ['1e2','1e3','1e4','1e5','1e6']
    values =  [1e2,1e3,1e4,1e5,1e6]
    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]

    gFlopE = []
    gFlopRK2 = []
    gFlopRK4 = []
    gFlopSE = []
    gFlopL =[]
    ar_data = [gFlopE, gFlopRK2, gFlopRK4, gFlopSE, gFlopL]


    plt.title("Orbit integration benchmark")


    for step in nsteps:
        bFile = open('%s/benchmark.txt'%(step))
        csvreader = csv.reader(bFile,delimiter=' ')

        header = []
        header = next(csvreader)

        rows = []
        count = 0
        for row in csvreader:
            rows.append(row)
            ar_data[count].append(row[1])
            count += 1

    for i in range(len(ar_data)):
        plt.loglog(values,ar_data[i],label=functions[i])

    plt.xlabel("# Steps")
    plt.ylabel("Performance [GFlop/s]")
    plt.legend()
    plt.savefig("plots/benchmark_plot.png",bbox_inches='tight')


main()




