Ce programme a été réalise dans le cadre de l'UE HAP924P - Modélisation Atomistique du parcours de Physique Numérique de la Faculté des Sciences

Il s'agit du travail commun de : AL-BASRI Ahmed, DATTIGNIE Enzo, EL-KOUIJAR Khadija.

# Utilisation du programme

Compiler le programme avec la commande : gcc -O2 main.c -o filename.out

Exécuter le programme avec la commande : ./filename.out 

Le programme accepte soit 
- aucun argument et les conditions de base sont appliquées
- un seul argument gérant la seed
- 3 arguments gérant respectivement : le nombre de pas jusqu'a un temps t_star, le mode d'intégration (0 de Euler, 1 pour verlet), et le temps t_star

Le programme crée un fichier nativement nommé "temp.txt" dans lequel il stocke 10000 séquence de mesures importantes notées dans le TP

Le programme exec.py permet la série de mesures menant aux graphiques de l'évolution de l'erreur en fonction de dt
