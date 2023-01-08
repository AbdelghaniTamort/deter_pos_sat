import numpy as np
import math as m
import textwrap

# ouvir le fichier de navigation et saut des lignes de commentaire
f = open("falk3440.21n", "r")
for line in f:
    if "END OF HEADER" in line:
        break


# extraction des données input exigées pour le calcul des coordonnées des satellites
list1 = []
for line in f:
    list1.append(line.strip())

f.close()


i = 0
list2 = []
b = len(list1)
while i < (len(list1)):
    list2.append(list1[i + 1])
    list2.append(list1[i + 2])
    list2.append(list1[i + 3])
    list2.append(list1[i + 4])
    list2.append(list1[i + 5])
    i = i + 8


list3 = []
i = 0
for i in list2:
    if len(i) == 75:
        b = []
        b.append(i[0:18])
        b.append(i[18:37])
        b.append(i[37:56])
        b.append(i[56:76])
        list3.append(b)

    if len(i) == 76:
        c = textwrap.wrap(i, 19)
        list3.append(c)


list4 = []
f = 0
while f < (len(list3)):
    list4.append(list3[f] + list3[f + 1] + list3[f + 2] +
                 list3[f + 3] + list3[f + 4])
    f = f + 5

list5 = []
for i in list3:
    n = [float(value.replace("D", "e")) for value in i]
    list5.append(n)


# choix du satellite selon son ordre dans le fichier Rinex de navigation
n_sat = int(
    input("veuillez entrer le numero d'ordre du satellite dans le fichier rinex: "))
list6 = []
min = (n_sat - 1) * 5
max = n_sat * 5
list6.append(list5[min:max])


# introduction du temps de mesure en format heures/minutes/secondes
t_h = float(input("temps de mesure en heure: "))
t_min = float(input("temps de mesure en minute: "))
t_s = float(input("temps de mesure en seconde: "))
t = t_h*3600+t_min*60+t_s

# attribution des données input à leurs variables respectives
t0e = list6[0][0][0]  # epoque origine
m0 = list6[0][0][3]  # Anomalie moyenne à l’époque origine
Delta_N = list6[0][0][2]  # Différence du mouvement moyen
a = list6[0][1][3]  # racine du demi grand axe
e = list6[0][1][1]  # excentricité
omega = list6[0][3][2]  # argument du périgée
# ascension droite du noeud ascendant à l'époque origine
omega_0 = list6[0][2][2]
cuc = list6[0][1][0]  # corrections cuc
cus = list6[0][1][2]  # corrections cus
cic = list6[0][2][1]  # corrections cic
crc = list6[0][3][1]  # corrections crc
crs = list6[0][0][1]  # corrections crs
cis = list6[0][2][3]  # corrections cis
i_point = list6[0][4][0]  # Taux de variation de l'angle d'inclinaison
i0 = list6[0][3][0]  # Inclinaison à  l'époque origine
omega_point = list6[0][3][3]  # vitesse de l'ascension droite

# Calculs
tk = t-t0e
if tk > 302400:
    tk = tk-604800
elif tk < -302400:
    tk = tk+604800
GM = 3.986004418*10e14  # en m**3s**-2
omega_e = 7.292115e-5
mk = m0+((np.sqrt(GM)/a**3)+Delta_N)*tk
Ek = mk
E = 0
while abs(Ek-E) > 1E-9:
    E = Ek
    Ek = Ek+e*m.sin(Ek)

# calcul de l'anomalie vraie
vk = m.atan((m.sqrt(1-e**2)*m.sin(Ek))/(m.cos(Ek)-e))
# calcul de l'argument de la latitude Uk
uk = omega + vk + cuc*m.cos(2*(omega+vk)) + cus*m.sin(2*(omega+vk))
rk = (a**2)*(1-e*m.cos(Ek))+crc*m.cos(2*(omega+vk))+crs * \
    m.sin(2*(omega+vk))  # calcul de la distance radiale rk
# calcul de l'inclinaison ik du plan orbital
ik = i0+i_point*tk+cic*m.cos(2*(omega+vk))+cis*m.sin(2*(omega+vk))
lambda_k = omega_0+(omega_point-omega_e)*tk-omega_e * \
    t0e  # longitude du noeud ascendant lambda_k

# definition des matrices de rotation
R3_lambda = np.array([[m.cos(-lambda_k), -m.sin(-lambda_k), 0],
                     [m.sin(-lambda_k), m.cos(-lambda_k), 0], [0, 0, 1]])
R3_uk = np.array([[m.cos(-uk), -m.sin(-uk), 0],
                 [m.sin(-uk), m.cos(-uk), 0], [0, 0, 1]])
R1_ik = np.array([[1, 0, 0], [0, m.cos(-ik), -m.sin(-ik)],
                 [0, m.sin(-ik), m.cos(-ik)]])

# Calcul des coordonnées dans le référentiel CTS
vect = np.array([[rk], [0], [0]])
prod = np.dot(np.dot(R3_lambda, R1_ik), R3_uk)
coord = np.dot(prod, vect)
print("les coordonnées du satellite N", n_sat, " sont: X=",
      coord[0], " et Y=", coord[1], " et Z=", coord[2])
