#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    resolutions = ["1e2","1e3","1e4","1e5","1e6"]
    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]

    for res in resolutions:
        plt.figure()
        plt.title(res)
        for f in functions:
            data = pd.read_csv("%s/%s.txt"%(res,f),sep=' ')

            plt.plot(data['X'],data['Y'],label=f)

        plt.axis("equal")
        plt.legend()
        plt.savefig("plots/%s_plot.png"%(res),bbox_inches='tight')
        #plt.show()

main()




