import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os


times = []
for i in range(2,5):
    for j in range(1,9):
        times.append(str(j)+"e"+str(i))
modes = ["0","1"] #0 = Euler, 1 = Verlet
t_star = "0.0160"

samples = int(1e5)
for mode in modes :
    ln_dE = []
    ln_dt = []
    dic = {}
    for time in times :
        # os.startfile("./test.out "+time+" "+mode+" "+t_star)
        subprocess.run(["./a.out",time,mode,t_star])
        with open ("./temp.txt") as f:
            compteur = -2
            dt = 0
            E = np.zeros(samples)
            t = np.zeros(samples)
            for Line in f :
                if (compteur==-2):
                    dt = Line.split("=")[1].strip()
                if compteur>=0 :
                    # if compteur < 10 :
                        # print(Line.split(";")[1].strip())
                    t[compteur] = Line.split(";")[0].strip()
                    E[compteur] = Line.split(";")[1].strip()
                compteur += 1
            E = E[0:compteur]
            t = t[0:compteur]
            E0 = E[0]
            for i in range(0,len(E)):
                E[i] = E[i] - E0
                if (t[i] == float(t_star)) :
                    ln_dE.append(np.log(abs(E[i])))
                    ln_dt.append(np.log(float(dt)))
                    print(f"dE = {E[i]}, dt = {dt}")
            t = t[:compteur]
            plt.plot(t,E)  
            dic[dt]=E
    plt.show()
    new_dE=[]
    new_dt=[]
    for i in range(0,len(ln_dt)):
        small = min(ln_dt)
        for i in range(0,len(ln_dt)) :
            if (ln_dt[i] == small) :
                new_dt.append(small)
                new_dE.append(ln_dE[i])
                ln_dt.pop(i)
                ln_dE.pop(i)
    if mode == "0" :
        temp1,temp2 = new_dt,new_dE
    else :
        plt.plot(temp1,temp2,"o",label="Euler")
        plt.plot(new_dt,new_dE,"o",label="Verlet")

a_euler = (temp2[12]-temp2[4])/(temp1[12]-temp1[4])
a_verlet = (new_dE[12]-new_dE[8])/(new_dt[12]-new_dt[8])
b_euler = temp2[4]-a_euler*temp1[4  ]
b_verlet = new_dE[8]-a_verlet*new_dt[8]

dt_copy = np.copy(new_dt)
droite_euler = a_euler*dt_copy+b_euler
droite_verlet = a_verlet*dt_copy+b_verlet


print(f"a_euler = {a_euler}, a_verlet = {a_verlet}")

plt.plot(dt_copy,droite_euler)
plt.plot(dt_copy,droite_verlet)


plt.legend()
plt.show()
