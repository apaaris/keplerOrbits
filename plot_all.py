#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]
    ec = '1e3'
    plt.figure()
    plt.title("Comparison of different methods")
    for f in functions:
        data = pd.read_csv("%s/%s.txt"%(ec,f),sep=' ')
        plt.plot(data['X'],data['Y'],label=f)

    plt.axis("equal")
    plt.legend()
    plt.savefig("plots/plot_all.png", bbox_inches='tight')
        #plt.show()

main()




