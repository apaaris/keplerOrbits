#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]
    #ec = '1e3'
    #ecs = ['1e2','1e3','1e4','1e5','1e6']
    ecs = ['1e4']
    for e in ecs:
        plt.figure()
        plt.title("Energy delta")

        for f in functions:
            data = pd.read_csv("%s/%s.txt_EM.txt"%(e,f),sep=' ')
            plt.plot(data['E'],label=f)

        plt.legend()
        plt.savefig("plots/plot_E.png", bbox_inches='tight')
    for e in ecs:
        plt.figure()
        plt.title("Momentum delta")

        for f in functions:
            data = pd.read_csv("%s/%s.txt_EM.txt"%(e,f),sep=' ')
            plt.plot(data['M'],label=f)

        plt.legend()
        #plt.show()
        plt.savefig("plots/plot_M.png", bbox_inches='tight')

main()




