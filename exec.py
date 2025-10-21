import subprocess
import matplotlib.pyplot as plt
import os


times = ["1e2","3e2"]
modes = ["0","1"]
t_star = "1"

samples = int(1e5)
for mode in modes :
    dic = {}
    for time in times :
        # os.startfile("./test.out "+time+" "+mode+" "+t_star)
        subprocess.run(["./test.out",time,mode,t_star])
        with open ("./temp.txt") as f:
            compteur = -2
            dt = 0
            E = 0
            t = 0
            for Line in f :
                if (compteur==-2):
                    dt = Line.split("=")[1].strip()
                if compteur>=0 :
                    t[compteur] = Line.split(";")[0].strip()
                    E[compteur] = Line.split(";")[1].strip()
                compteur += 1
            E = E[:compteur]
            t = t[:compteur]
            plt.plot(t,E)
            plt.show()
            dic[dt]=E
     
    