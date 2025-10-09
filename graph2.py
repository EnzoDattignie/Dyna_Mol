import os
import numpy as np
import matplotlib.pyplot as plt

file = "./temp.txt"
samples = 
n_star = 1


t = np.zeros(samples)
dic = {}

with open (direc+"/"+files[0]) as f:
    compteur = -2
    dt = 0
    E = np.zeros(samples)
    for Line in f :
        if (compteur==-2):
            dt = Line.split("=")[1].strip()
        if compteur>=0 :
            t[compteur] = Line.split(";")[0].strip()
            E[compteur] = Line.split(";")[1].strip()
        compteur += 1
    dic[dt]=E
    
for file in files[1:] :
    with open (direc+"/"+file) as f:
        compteur = -2
        dt = 0
        E = np.zeros(samples)
        for Line in f :
            if (compteur==-2):
                dt = Line.split("=")[1].strip()
            if compteur>=0 :
                E[compteur] = Line.split(";")[1].strip()
            compteur += 1
            dic[dt]=E

Liste_dt=np.zeros(len(files))
Liste_dE=np.zeros(len(files))

compteur = 0
for elmt in dic :
    # plt.plot(t,dic[elmt])
    # plt.ylim(5,7)
    Liste_dt[compteur]=elmt
    Liste_dE[compteur]=abs(dic[elmt][n_star]-dic[elmt][0])
    compteur += 1

# plt.show()
for i in range (0,len(Liste_dt)-1):
    for j in range (i+1,len(Liste_dt)) :
        if Liste_dt[i]>Liste_dt[j] :
            temp = Liste_dt[j]
            temp2 = Liste_dE[j]
            Liste_dt[j] = Liste_dt[i]
            Liste_dE[j] = Liste_dE[i]
            Liste_dt[i] = temp
            Liste_dE[i] = temp2

ln_dt1 = np.log10(Liste_dt)
ln_dE1 = np.log10(Liste_dE)

direc = "./Euler2"
files = os.listdir(direc)
samples = 51
n_star = 20


t = np.zeros(samples)
dic = {}

with open (direc+"/"+files[0]) as f:
    compteur = -2
    dt = 0
    E = np.zeros(samples)
    for Line in f :
        if (compteur==-2):
            dt = Line.split("=")[1].strip()
        if compteur>=0 :
            t[compteur] = Line.split(";")[0].strip()
            E[compteur] = Line.split(";")[1].strip()
        compteur += 1
    dic[dt]=E
    
for file in files[1:] :
    with open (direc+"/"+file) as f:
        compteur = -2
        dt = 0
        E = np.zeros(samples)
        for Line in f :
            if (compteur==-2):
                dt = Line.split("=")[1].strip()
            if compteur>=0 :
                E[compteur] = Line.split(";")[1].strip()
            compteur += 1
            dic[dt]=E

Liste_dt=np.zeros(len(files))
Liste_dE=np.zeros(len(files))

compteur = 0
for elmt in dic :
    # plt.plot(t,dic[elmt])
    # plt.ylim(5,7)
    Liste_dt[compteur]=elmt
    Liste_dE[compteur]=abs(dic[elmt][n_star]-dic[elmt][0])
    compteur += 1

# plt.show()
for i in range (0,len(Liste_dt)-1):
    for j in range (i+1,len(Liste_dt)) :
        if Liste_dt[i]>Liste_dt[j] :
            temp = Liste_dt[j]
            temp2 = Liste_dE[j]
            Liste_dt[j] = Liste_dt[i]
            Liste_dE[j] = Liste_dE[i]
            Liste_dt[i] = temp
            Liste_dE[i] = temp2

ln_dt = np.log10(Liste_dt)
ln_dE = np.log10(Liste_dE)

# print(f"Fit de la forme {alpha}*x + {b}")
data1, = plt.plot(ln_dt1,ln_dE1,'o',label="Verlet")
data2, = plt.plot(ln_dt,ln_dE,'o',label="Euler")
plt.title("Graphique représentant la différence dE en fonction de dt")
plt.xlabel("ln(dt)")
plt.ylabel("ln(dE)")
plt.legend()


# x = np.linspace(-14,0,50)
# alpha = (ln_dE[int(len(ln_dE1)/2)+2]-ln_dE[int(len(ln_dE)/2)-2])/(ln_dt[int(len(ln_dE)/2)+2]-ln_dt[int(len(ln_dE)/2)-2])
# b = ln_dE[int(len(ln_dE)/2)]-alpha*ln_dt[int(len(ln_dE)/2)]
# y = alpha*x+b
# plt.plot(x,y)


plt.show()

# for file in files :
#     with open (file) as f :
#         for lines in f: