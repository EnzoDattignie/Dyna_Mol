import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os


times = ["1","5","10","50","100"]
modes = ["0","1"]
t_star = "0.06"

samples = int(1e5)
for mode in modes :
    dic = {}
    for time in times :
        # os.startfile("./test.out "+time+" "+mode+" "+t_star)
        subprocess.run(["./test.out",time,mode,t_star])
        with open ("./temp.txt") as f:
            compteur = -2
            dt = 0
            E = np.zeros(samples)
            t = np.zeros(samples)
            for Line in f :
                if (compteur==-2):
                    dt = Line.split("=")[1].strip()
                if compteur>=0 :
                    if compteur < 10 :
                        print(Line.split(";")[1].strip())
                    t[compteur] = Line.split(";")[0].strip()
                    E[compteur] = Line.split(";")[1].strip()
                compteur += 1
            E = E[0:compteur]
            t = t[0:compteur]
            E0 = E[0]
            for i in range(0,len(E)):
                E[i] = E[i] - E0
            t = t[:compteur]
            plt.plot(t,E)
            plt.show()
            dic[dt]=E
     