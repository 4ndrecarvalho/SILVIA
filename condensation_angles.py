#! /usr/bin/env python3
#_________________________________________________________________________________
# andre.dc@ua.pt
# MolModel@CICECO Universidade Aveiro
#
# This script calculate the mol fraction of silica species in the simulation, and generate the graphs for the (mole fraction / time) and (mole fraction / degree of condensation)
#
# USAGE: REP=<insert_Erep> && gmx trjconv -nobackup -s $REP.tpr -f $REP.xtc -o d.xtc -n ../index.ndx  && python3 ../condensation_legacy.py d.xtc <insert_distance> no $REP
#
# <insert_Erep> is the value for VS repulsion energy, the index.ndx contains only Si atoms and <insert_distance> is the radial distance to count a bond, 5.2  for SQSi, 5.8 for QSi.
#__________________________________________________________________________________


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis
import sys
import tqdm

input = sys.argv[1]
cutoff = np.array(float(sys.argv[2]))
av = sys.argv[3]
fname = sys.argv[4] if len(sys.argv) == 5 else ""
system = MDAnalysis.Universe(input)
n_dummies = 4  # número de dummies por partícula
n_atoms = (len(system.atoms)-1000)/n_dummies    # (len(total_particles)- silica particles)
nb_dummies = []
bonds = []
dij = []
overall = []
d_nb = []
si_nb = []
time = []

def dist_PBC(x0, x1, dimensions):
    # https://stackoverflow.com/questions/11108869/optimizing-python-distance-calculation-while-accounting-for-periodic-boundary-co
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def mov_avg(data, time, window):
    avg = [np.mean(data[i:i+window]) for i in range(len(data)-window)]
    t_avg = time[window//2:-window//2]
    return t_avg, avg


def calc_angle(si1, d1, si2, d2):
    v1 = d1 - si1
    v1 /= np.linalg.norm(v1)
    v2 = d2 - si2
    v2 /= np.linalg.norm(v2)
    d1 = np.sqrt((v1[0])**2 + (v1[1])**2 + (v1[2])**2)
    d2 = np.sqrt((v2[0])**2 + (v2[1])**2 + (v2[2])**2)
    return np.arccos(np.dot(v1,v2)/(d1*d2))


with open("condensation_Si-D-D-Si_angles_"+fname+".xvg", "a+") as flog:
    print("# Andre Carvalho\n# MolModel Group\n# Universidade de Aveiro\n# andre.dc@ua.pt", file=flog)                                                          
    print('@    title "Silica condensation Si-D D-Si angles"\n@    xaxis  label "time"\n@    yaxis  label "mean angle ( º )"\n@TYPE xy\n@ view 0.15, 0.15, 0.75, 0.85\n', file=flog)          
    with open("condensation_dummy_"+fname+".xvg", "w+") as fout:
        print("# Andre Carvalho\n# MolModel Group\n# Universidade de Aveiro\n# andre.dc@ua.pt", file=fout)
        print('@    title "Silica condensation"\n@    xaxis  label "time"\n@    yaxis  label "qi"\n@TYPE xy\n@ view 0.15, 0.15, 0.75, 0.85\n@ legend on\n@ legend box on\n@ legend loctype view\n@ legend 0.78, 0.8\n@ legend length 2\n@ s0 legend "Q0"\n@ s1 legend "Q1"\n@ s2 legend "Q2"\n@ s3 legend "Q3"\n@ s4 legend "Q4"\n@ s5 legend "Qn"\n@ s6 legend "C"\n@ s0 line color "black"\n@ s1 line color "red"\n@ s2 line color "blue"\n@ s3 line color "magenta"\n@ s4 line color "brown"\n@ s5 line color "yellow"\n@ s6 line color "green"', file=fout)
        for t, ts in zip(tqdm.tqdm(range(system.trajectory.n_frames)), system.trajectory):
            box = system.dimensions
            time.append(ts.time)
            coords = system.atoms.positions
            si_coords = coords[[i for i in range(len(coords)) if i %5 == 0]]
            d_coords = coords[[i for i in range(len(coords)) if i %5 != 0]]
            # criar matriz de distâncias de cada dummy a todos os dummies
            [dij.append(dist_PBC(d_coords, i, box[0:3])) for i in d_coords]
            # cria matriz de boleanos, True se o dummy estiver abaixo do cutoff
            near_dummies = (dij < cutoff)
            # remove a interacção com os dummies intramoleculares
            for i in range(0, len(near_dummies), n_dummies):
                near_dummies[i:i+n_dummies, i:i+n_dummies] = False
            # criar lista de listas com identificação dos vizinhos de cada dummie
            nb_dummies = [np.where(line) for line in near_dummies]
    
            [d_nb.append((d[0])) for d in nb_dummies]
            angle = 0
            ang_count = 0
            dist_si_si = 0
            for i,d in enumerate(d_nb):
                if len(d) == 0:
                    continue
                elif len(d) >= 1:
                    for j in d:
                        angle += calc_angle(si_coords[int(np.floor(i/4))], d_coords[i], si_coords[int(np.floor(j/4))], d_coords[int(j)])
                        dist_si_si += dist_PBC(si_coords[int(np.floor(i/4))],si_coords[int(np.floor(j/4))], box[0:3])
                        ang_count += 1

            print("{} {:.3f}  {:.3f}".format(ts.time, angle/ang_count*180/np.pi, dist_si_si/ang_count), file=flog)        

            # Colapsar vizinhos de %n_dummies para tornar a lista de listas de vizinhos de cada dummy em lista de listas de vizinhos de cada partícula
            for _ in range(0, len(nb_dummies), n_dummies):
                # contabilizar vizinhança entre particula e dummies ligados
                # d_nb.append([dummies_near for i in [_] for particle in nb_dummies[i:i+n_dummies] for dummies_near in np.array(particle).flatten()]) 

                # contabilizar vizinhança de acordo com a que particulas os dummies ligados pertencem
                si_nb.append([np.floor(dummies_near/n_dummies) for i in [_] for particle in nb_dummies[i:i+n_dummies] for dummies_near in np.array(particle).flatten()]) 

            # a primeira frame não tem informação de estado anterior
            if ts.frame == 0:
                nb_old = si_nb.copy()
            else:
                # comparar lista de vizinhos actual com a anterior
                for old, new in zip(nb_old, si_nb):
                    # se vizinho é o mesmo então está ligado
                    bonds.append([i for i in old for j in new if i == j])
                nb_old = si_nb.copy()
                # cria lista com o número de ligações de cada átomo
                overall.append([len(x) for x in bonds])
                Q0 = Q1 = Q2 = Q3 = Q4 = Qn = 0
                for i in [len(x) for x in bonds]:
                    if i == 0:
                        Q0 += 1
                    elif i == 1:
                        Q1 += 1
                    elif i == 2:
                        Q2 += 1
                    elif i == 3:
                        Q3 += 1
                    elif i == 4:
                        Q4 += 1
                    elif i >= 5:
                        Qn += 1
                c = 1/4*((Q1/n_atoms+2*(Q2/n_atoms)+3*(Q3/n_atoms)+4*(Q4/n_atoms)))
                print("{} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}".format(
                    time[t], Q0/n_atoms, Q1/n_atoms, Q2/n_atoms, Q3/n_atoms, Q4/n_atoms, Qn/n_atoms, c), file=fout)
                Q0 = Q1 = Q2 = Q3 = Q4 = Qn = 0
                bonds.clear()
            nb_dummies.clear()
            si_nb.clear()
            d_nb.clear()
            dij.clear()
# ### MAKE GRAPH #####
q0 = []
q1 = []
q2 = []
q3 = []
q4 = []
qn = []
C = []
Q0 = Q1 = Q2 = Q3 = Q4 = Qn = 0
for frame in overall:
    for i in frame:
        if i == 0:
            Q0 += 1
        elif i == 1:
            Q1 += 1
        elif i == 2:
            Q2 += 1
        elif i == 3:
            Q3 += 1
        elif i == 4:
            Q4 += 1
        elif i >= 5:
            Qn += 1
    c = 1/4*((Q1/n_atoms)+2*(Q2/n_atoms)+3*(Q3/n_atoms)+4*(Q4/n_atoms))
    q0.append(Q0/n_atoms)
    q1.append(Q1/n_atoms)
    q2.append(Q2/n_atoms)
    q3.append(Q3/n_atoms)
    q4.append(Q4/n_atoms)
    qn.append(Qn/n_atoms)
    C.append(c)
    Q0 = Q1 = Q2 = Q3 = Q4 = Qn = 0
with open("condensation_dummy_"+fname+".log", "a+") as fout:
    print("\nQ0: {:3.3f}% | Q1: {:3.3f}% | Q2: {:3.3f}% | Q3: {:3.3f}% | Q4: {:3.3f}% \nQ5+: {:3.3f}% | C:  {:3.3f}% | Max CN: {}  | Mean CN: {:3.3f}".format(np.mean(q0)
                                                                                                                                                              * 100, np.mean(q1)*100, np.mean(q2)*100, np.mean(q3)*100, np.mean(q4)*100, np.mean(qn)*100, np.mean(C)*100, np.max(overall), np.mean(overall)), file=fout)
if system.trajectory.n_frames-1 >= 1000 and av == True:  # usar médias móveis para trajectórias maiores
    window = (system.trajectory.n_frames-1)//100
    _, q0 = mov_avg(q0, time[1:], window)
    _, q1 = mov_avg(q1, time[1:], window)
    _, q2 = mov_avg(q2, time[1:], window)
    _, q3 = mov_avg(q3, time[1:], window)
    _, q4 = mov_avg(q4, time[1:], window)
    _, qn = mov_avg(qn, time[1:], window)
    x, C = mov_avg(C, time[1:], window)
else:
    x = time[1:]
l = np.log10(x)

plt.rcParams["font.family"] = "Times New Roman"
fig = plt.figure()
ax = plt.subplot(111)
plt.plot(x, q0, c='black', label=r"Q$_0$")  # , marker="+")
plt.plot(x, q1, c='#EE1D23', label=r"Q$_1$")  # , marker="s")
plt.plot(x, q2, c='#274FA2', label=r"Q$_2$")  # , marker="D")
plt.plot(x, q3, c='#8C4198', label=r"Q$_3$")  # , marker="^")
plt.plot(x, q4, c='#6C0C40', label=r"Q$_4$")  # , marker="o")
plt.plot(x, C, c='#008D48', label=r"C")  # , marker="v")
plt.plot(x, qn, c='yellow', label=r"Q$_{5^+}$")
plt.title(fname)
plt.xlabel("time")
plt.ylabel("qi")
plt.ylim(0, 1)
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("condensation_dummy_"+fname+".svg",
            bbox_extra_artists=[lgd], bbox_inches='tight')

fig = plt.figure()
ax = plt.subplot(111)
plt.plot(l, q0, c='black', label=r"Q$_0$")  # , marker="+")
plt.plot(l, q1, c='#EE1D23', label=r"Q$_1$")  # , marker="s")
plt.plot(l, q2, c='#274FA2', label=r"Q$_2$")  # , marker="D")
plt.plot(l, q3, c='#8C4198', label=r"Q$_3$")  # , marker="^")
plt.plot(l, q4, c='#6C0C40', label=r"Q$_4$")  # , marker="o")
plt.plot(l, C, c='#008D48', label=r"C")  # , marker="v")
plt.plot(l, qn, c='yellow', label=r"Q$_{5^+}$")
plt.title(fname)
plt.xlabel("time")
plt.ylabel("qi")
plt.ylim(0, 1)
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("condensation_dummy_"+fname+"log.svg",
            bbox_extra_artists=[lgd], bbox_inches='tight')

fig = plt.figure()
ax = plt.subplot(111)
plt.plot(C, q0, c='black', label=r"Q$_0$")  # , marker="+")
plt.plot(C, q1, c='#EE1D23', label=r"Q$_1$")  # , marker="s")
plt.plot(C, q2, c='#274FA2', label=r"Q$_2$")  # , marker="D")
plt.plot(C, q3, c='#8C4198', label=r"Q$_3$")  # , marker="^")
plt.plot(C, q4, c='#6C0C40', label=r"Q$_4$")  # , marker="o")
plt.plot(C, qn, c='yellow', label=r"Q$_{5^+}$")
plt.title(fname)
plt.xlabel("Degree of condensation")
plt.ylabel("qi")
plt.ylim(0, 1)
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("condensation_dummy_"+fname+"C.svg",
            bbox_extra_artists=[lgd], bbox_inches='tight')
