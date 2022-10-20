#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    resolutions = ["1e4","1e5","1e6","1e7"]
    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]

    for f in functions:
        plt.figure()
        plt.title(f)
        for res in resolutions:
            data = pd.read_csv("%s/%s.txt"%(res,f),sep=' ')

            plt.plot(data['X'],data['Y'],label=res)

        plt.axis("equal")
        plt.legend()
        plt.savefig("plots/%s_plot.png"%(f),bbox_inches='tight')
        #plt.show()

main()




