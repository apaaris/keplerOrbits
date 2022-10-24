#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():

    e = ['ec1','ec2','ec3','ec4','ec5','ec6','ec7','ec8','ec9']
    functions = ["eEuler", "rk2","rk4","semiI","leapfrog"]

    for ec in e:
        plt.figure()
        plt.title(ec)
        for f in functions:
            data = pd.read_csv("%s/%s.txt"%(ec,f),sep=' ')

            plt.plot(data['X'],data['Y'],label=f)

        plt.axis("equal")
        plt.legend()
        plt.savefig("plots/%s_plot.png"%(ec),bbox_inches='tight')
        #plt.show()

main()




