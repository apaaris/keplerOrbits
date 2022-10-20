#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    e = ['ec1','ec2','ec3','ec4','ec5','ec6','ec7','ec8','ec9']
    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]

    for f in functions:
        plt.figure()
        plt.title(f)
        for ec in e:
            data = pd.read_csv("%s/%s.txt"%(ec,f),sep=' ')

            plt.plot(data['X'],data['Y'],label=ec)

        plt.axis("equal")
        plt.legend()
        plt.savefig("plots/%s_plot_ec.png"%(f),bbox_inches='tight')
        #plt.show()

main()




