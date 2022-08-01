import math
from numpy import array, zeros, fabs, linalg
import numpy as np
import time
import random
from numpy import linalg as LA
print("----Pivot de Gauss Totale----")
print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")
print("vous devez choisir:")
print("1) votre propre matrice")
print("2) une matrice aléatoire ")
choix=int(input("choisir 1 / 2  "))
while(not(choix in range(1,3))):
    choix=int(input("choisir 1 / 2   "))


if (choix==1) :
    R = int(input("Entrer la taille de matrice A:"))
    while(not(R in range(2,100))):
        R = int(input("Entrer la taille de matrice A:"))

    print("donner la Matrice A")
    entries = list(map(float, input().split()))

    # For printing the matrix
    a = np.array(entries).reshape(R, R)
    print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

    print("A= ")
    print(a)
    print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

    print("donner le vecteur B")
    entries1 = list(map(float, input().split()))

    # For printing the matrix
    b = np.array(entries1).reshape(R, 1)
    print("b= ")
    print(b)

else:
    a = np.random.uniform(0, 20, size=(5, 5))

    b = np.random.uniform(0, 20, size=(5, 1))

    print(a)
    print(b)
condd=LA.cond(a)
n = len(b)

solut = np.linalg.solve(a,b)    #on va utiliser numpy pour vérifier qu'on a trouver la solution correcte
print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

print("----la solution avec numpy----")
print(solut)

x=zeros(n,float)







def zero(a,b,n):     # la fonction zero permet de mettre on zero les valeur qui sont négligable c à d qui sont trés petits
    eps=1^(-4)
    for i in range (0,n):
        for j in range (0,n):
            if eps>(fabs(a[i][j])):
                a[i][j] = 0

        if eps>(fabs(b[i])):
            b[i]=0




start = time.time()     #pour calculer le temps d'execution on fait commencer un timer au début de traitement



pivot_sol=zeros(n,int)

for i in range (0,n):  #vecteur de pivotation des solutions
    pivot_sol[i]=i




for k in range(0,(n-1)):
    ref=0
    for i in range(k,n):

        for j in range(k,n):
            if ref<fabs(a[i][j]):                         #chercher le pivot max dans la matrice A
                ref=fabs(a[i][j])
                ligne=i
                colonne=j                                  #chaque fois en trouve une valeur supérieur a la valeur précedante
                                                           #on garde sa ligne et sa colonne


    for j in range (k,n):       #une fois qu'on trouve un pivot max on fait la permutation des ligne pour la matrice A
        temp=a[k][j]
        a[k][j] = a[ligne][j]
        a[ligne][j]=temp


    b[k][0],b[ligne][0]=b[ligne][0],b[k][0]


    for i in range(0,n):  #les ligne sont permutés mainteneant on fait la permutation des colonnes
        temp=a[i][k]
        a[i][k]=a[i][colonne]
        a[i][colonne]=temp

    temps=pivot_sol[k]
    pivot_sol[k]=pivot_sol[colonne]
    pivot_sol[colonne]=temps

    if (a[k][k] == 0):   #une fois que le pivot max est null alors on ne peux pas faire la méthode de pivot total
        print(" Un pivot nul ! => methode de Gauss pivot total non applicable")


    for i in range(k+1,n):
        p = (a[i][k]/a[k][k])  #maintenant on fait la méthode de gauss exemple L2 --> L2-a11/a22*L1
        for j in range(k,n):
            a[i][j]=a[i][j]-p*a[k][j]


        b[i]=b[i]-p*b[k]
print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

print("----la matrice A aprés la modification----")
print(a)
print("----le vecteur B aprés la modification----")
print(b)
for i in reversed(range(n)):
    s=0
    for j in range (i+1,n):
        s=s+a[i][j]*b[j]      #ce boucle nous permet de determiner les solutions x,y,z,..
    b[i]=(b[i]-s)/a[i][i]

for i in range(0,n):    #l'association de chaque valeurs a sa soltuion exemple :
                        # X=[z,x,y] ========> X=[x,y,z]
    x[pivot_sol[i]]=b[i]




zero(a,b,n)
end = time.time()
print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

print("le temps de calcul =", end - start)
print("cond(A)=",condd)
print("<*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*><*>")

print("---la solution du systéme avec Pivot de Gauss Total ---")

print(x)

del a
del b
del x











