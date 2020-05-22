#! /usr/bin/env python3
import numpy as np
import sys
import os

input=sys.argv[1]
# directory = "/home/acarvalho/Desktop/silvia/md/stella8VS_new_topology_Dall/"

fname=sys.argv[1]
#Valores experimentais de Devreux et al.
#fracção molar (q) e grau de condensação (c) --> (q, c)
d_max_c       = .9
d_max_q0      = np.array([1.0, .00])
d_max_q1      = np.array([.58, .25])
d_max_q2      = np.array([.58, .52])
d_max_q3      = np.array([.60, .77])
d_max_q4      = np.array([.47, .90])
d_cross_q0_q1 = np.array([.45, .17])
d_cross_q1_q2 = np.array([.45, .41])
d_cross_q2_q3 = np.array([.44, .64])
d_cross_q3_q4 = np.array([.50, .87])
d_cross_q0_q2 = np.array([.20, .26])
d_cross_q1_q3 = np.array([.22, .32])
d_cross_q2_q4 = np.array([.22, .75])
d_cross_q0_q3 = np.array([.06, .36])
d_cross_q1_q4 = np.array([.07, .62])
d_cross_q0_q4 = np.array([.01, .46])

def distance(i, j):
    return np.sqrt((i[0] - j[0])**2 + (i[1] - j[1])**2)

def crossing_points_list(a, b, c, cutoff, name):
    """
    Avalia os pontos de cruzamento entre listas
    cutoffs foram aplicados para remover cruzamentos
    inuteis, como os das formaçoes iniciais de q3 e q4
    ou os cruzamentos finais entre especies q0 e q1.
    """
    q_value = c_index = c_cross_a_b = 0
    a=np.array(a)
    b=np.array(b)
    a[a<cutoff]=0.001
    b[b<cutoff]=0.0001
    diff=(b-a)
    indexes = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]+1 #indexes where sign changes
    if len(indexes)<1:
        q_value = 0 
        c_cross_a_b = 0
        # print("No intersection in",name)
    else:
        c_index = int(np.round(np.mean(indexes)))
        c_cross_a_b = c[c_index]
        for i in indexes:
            q_value += (a[i] + b[i] + a[i-1] + b[i-1])/4
        q_value = q_value/len(indexes)
    # print("{:s} >> q_value: {:.2f}%, c_cross_a_b: {:.2f}%, indexes: {:d}, c_index: {:d}".format(name,q_value*100,c_cross_a_b*100,len(indexes),c_index))
    return [q_value, c_cross_a_b]

def xvg_reader(fname):
    """
    Leitor de xvg de condensação obtem informação
    sobre os máximos de condensação de cada espécio de Qn
    e os pontos de cruzamento entre espécies.
    """
    t  = []
    q0 = []
    q1 = []
    q2 = []
    q3 = []
    q4 = []
    qx = []
    c  = []
    with open(fname, 'r') as fin:
        for line in fin:
            if line.startswith(('#','@')):
                pass
            else:
                t.append(float(line.split()[0]))
                q0.append(float(line.split()[1]))
                q1.append(float(line.split()[2]))
                q2.append(float(line.split()[3]))
                q3.append(float(line.split()[4]))
                q4.append(float(line.split()[5]))
                qx.append(float(line.split()[6]))
                c.append(float(line.split()[7]))
        max_c = np.max(c)
        max_q0 = (np.max(q0), c[q0.index(np.max(q0))])
        max_q1 = (np.max(q1), c[q1.index(np.max(q1))]) 
        max_q2 = (np.max(q2), c[q2.index(np.max(q2))]) 
        max_q3 = (np.max(q3), c[q3.index(np.max(q3))]) 
        max_q4 = (np.max(q4), c[q4.index(np.max(q4))]) 
        max_qx = (np.max(qx), c[qx.index(np.max(qx))]) 
        cross_q0_q1 = np.array(crossing_points_list(q0, q1, c, .00, "q0_q1"))
        cross_q1_q2 = np.array(crossing_points_list(q1, q2, c, .00, "q1_q2"))
        cross_q2_q3 = np.array(crossing_points_list(q2, q3, c, .00, "q2_q3"))
        cross_q3_q4 = np.array(crossing_points_list(q3, q4, c, .04, "q3_q4")) #filtrar formações iniciais
        cross_q0_q2 = np.array(crossing_points_list(q0, q2, c, .00, "q0_q2"))
        cross_q1_q3 = np.array(crossing_points_list(q1, q3, c, .00, "q1_q3"))
        cross_q2_q4 = np.array(crossing_points_list(q2, q4, c, .00, "q2_q4"))
        cross_q0_q3 = np.array(crossing_points_list(q0, q3, c, .00, "q0_q3"))
        cross_q1_q4 = np.array(crossing_points_list(q1, q4, c, .00, "q1_q4"))
        cross_q0_q4 = np.array(crossing_points_list(q0, q4, c, .00, "q0_q4"))
    return max_c, max_q0, max_q1, max_q2, max_q3, max_q4, max_qx, cross_q0_q1, cross_q1_q2, cross_q2_q3, cross_q3_q4, cross_q0_q2, cross_q1_q3, cross_q2_q4, cross_q0_q3, cross_q1_q4, cross_q0_q4

def fitness(fname, qx_penalty=5.):
    """ 
    Descriptor de aptidão é baseado nos valores de
    condensação de cada espécie em função do grau
    de condensação(c) e da fracção molar (q).
    Máximos de cada espécie, pontos de cruzamento
    entre espécies.
    O cálculo é feito de acordo com a distancia
    euclideanea de cada ponto obtido ao de Devreux,
    com a adição de uma penalização por espécies
    com mais do que 4 ligações.
    """
    max_c, max_q0, max_q1, max_q2, max_q3, max_q4, max_qx, cross_q0_q1, cross_q1_q2, cross_q2_q3, cross_q3_q4, cross_q0_q2, cross_q1_q3, cross_q2_q4, cross_q0_q3, cross_q1_q4, cross_q0_q4 = xvg_reader(fname)    
    diff_max = distance(max_q0,d_max_q0) + distance(max_q1,d_max_q1) + distance(max_q2,d_max_q2) + distance(max_q3,d_max_q3) + distance(max_q4,d_max_q4)
    diff_cross = distance(cross_q0_q1,d_cross_q0_q1) + distance(cross_q1_q2,d_cross_q1_q2) + distance(cross_q2_q3,d_cross_q2_q3) + distance(cross_q3_q4,d_cross_q3_q4) + distance(cross_q0_q2,d_cross_q0_q2) + distance(cross_q1_q3,d_cross_q1_q3) + distance(cross_q2_q4,d_cross_q2_q4) + distance(cross_q0_q3,d_cross_q0_q3) + distance(cross_q1_q4,d_cross_q1_q4) + distance(cross_q0_q4,d_cross_q0_q4) 
    f = diff_max + diff_cross + qx_penalty*max_qx[0]
    print(f,"  >>>  ",fname)
    return f

class Model(object):
    constraints = []
    def __init__(self, chromosome, fitness, path="", generation=0):
        self.chromosome = chromosome
        self.fitness = fitness
        self.path = path
        self.generation = generation

    @staticmethod
    def from_xvg(fpath):
        parms = fpath.split("/")
        e_rep = float(parms[-1].split("condensation_")[1].split("rep_")[0])
        e_atract = float(parms[-2].split("_")[1])
        sigma = float(parms[-2].split("_")[0])
        d_vs = float(parms[-1].split("condensation_")[1].split("rep_")[-1].split("dvs.xvg")[0])
        # dt = float(path[-3].split("_")[1])   #only for final dirs.
        chromosome = [sigma, d_vs, e_atract, e_rep]
        fit = fitness(fpath)
        return Model(chromosome, fit, path=fpath)
    
    @staticmethod
    def from_log(fpath):
        population = []
        with open(fpath, 'r') as fin:
            for line in fin:
                try:
                    if ">>>" in line:
                        parm = line.split(">>>")
                        fitness = float(parm[0])
                        # try: # try to recalculate fitness
                        #     fitness = fitness(str(parm[1][:-1]))
                        # except: # get fitness from log file
                        fitness = float(parm[0])
                        print(fitness)
                        parms = parm[1].split("/")
                        e_rep = float(parms[-1].split("condensation_")[1].split("rep_")[0])
                        e_atract = float(parms[-2].split("_")[1])
                        sigma = float(parms[-2].split("_")[0])
                        d_vs = float(parms[-1].split("condensation_")[1].split("rep_")[-1].split("dvs.xvg")[0])
                        # dt = float(path[-3].split("_")[1])   #only for final dirs.
                        chromosome = [sigma, d_vs, e_atract, e_rep]
                        population.append(Model(chromosome, fitness, path=fpath))
                except:
                    pass
        return population

def run_md(child, directory="/home/acarvalho/Desktop/silvia/simplex/"):
    path = directory+"dVS_"+str(child.chromosome[1])+"/"+str(child.chromosome[0])+"_"+str(child.chromosome[2])+"/"
    # os.makedirs(path, exist_ok=True)
    os.system("bash "+directory+"run "+str(child.chromosome[0])+" "+str(child.chromosome[2])+" "+str(child.chromosome[3])+" "+str(child.chromosome[1]))
    os.system("gmx -nocopyright -nobackup trjconv -f "+path+str(child.chromosome[3])+".xtc -s "+path+str(child.chromosome[3])+".tpr -o "+path+"d.xtc -pbc mol -skip 2 -n "+directory+"index.ndx")
    os.system("python3 "+directory+"condensation.py "+path+"d.xtc 3.2 no "+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs") #Check if Qsi=5.8, SQsi=5.2, XQsi=3.2
    os.system("mv conde*rep* "+path)
    try:
        child.fitness = fitness(path+"condensation_"+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs.xvg")
    except:
        child.fitness = 999
    os.system("echo Fitness: "+str(child.fitness)+" >> "+path+"gene_dVS_"+str(child.chromosome[1])+"_"+str(child.chromosome[0])+"_"+str(child.chromosome[2])+"_"+str(child.chromosome[3])+".log")
    child.path = path+"condensation_"+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs.xvg"
    return child.chromosome, child.fitness, child.path

def simplex_nelder_mead(population, n_iter=10, n_vars=0, d_r=1, d_e=2, d_oc=.5, d_ic=-.5, d_s=.5):
    population.sort(key=lambda i: i.fitness, reverse=False)
    with open("pop_history.log", "a") as fout:
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", file=fout)
        print("Initial Population of "+str(input)+": ", file=fout)
        for model in population:
            print(str(model.chromosome)+"\t>->\t"+str(model.fitness), file=fout)
    simplex_history = []
    if n_vars == 0:
    # automatizar o n_vars através da comparação de valores entre modelos.
        for i in population:
            for j in population:
                diff = len(np.where((np.array(i.chromosome)-np.array(j.chromosome))!=0)[0])
                if diff > n_vars:
                    n_vars = diff
        n_vars += 1
    with open("pop_history.log", "a") as fout:
        print("\nNumber of simplex points: "+str(n_vars)+"\nNumber of iterations: "+str(n_iter), file=fout)
    for iter in range(0, n_iter):
        # Order
        population.sort(key=lambda i: i.fitness, reverse=False)
        population = population[:n_vars]
        with open("pop_history.log", "a") as fout:
            print("-------------------------------------------------------------------------------", file=fout)
            print("Iteration: "+str(iter), file=fout)
            print("-------------------------------------------------------------------------------", file=fout)
            print("Simplex points:", file=fout)
            print([(i.chromosome, i.fitness) for i in population], file=fout)
        centroid = [.0, .0, .0, .0]
        for i in population[:-1]:
            centroid = np.add(centroid, i.chromosome)
        centroid = centroid/len(population[:-1])

        # Reflect
        reflect = centroid + d_r * (centroid - population[-1].chromosome)
        _, r_fitness, r_path = run_md(Model(reflect, 0))
        simplex_history.append((r_fitness, r_path))
        with open("pop_history.log", "a") as fout:
            print("New reflection ("+str(reflect)+") fitness: "+str(r_fitness), file=fout)
        if r_fitness >= population[0].fitness and r_fitness < population[-1].fitness:
            population.pop(-1)
            population.append(Model(reflect, r_fitness))
            with open("pop_history.log", "a") as fout:
                print("Reflected ("+str(reflect)+"), fitness: "+str(r_fitness), file=fout)

        # Expand
        elif r_fitness < population[0].fitness:
            expand = centroid + d_e * (centroid - population[-1].chromosome)
            _, e_fitness, e_path = run_md(Model(expand, 0))
            if e_fitness <= r_fitness:
                population.pop(-1)
                population.append(Model(expand, e_fitness))
                simplex_history.append((e_fitness, e_path))
                with open("pop_history.log", "a") as fout:
                    print("Expanded ("+str(expand)+")", file=fout)
            else:
                population.pop(-1)
                population.append(Model(reflect, r_fitness))
                with open("pop_history.log", "a") as fout:
                    print("Reflected ("+str(reflect)+")", file=fout)

        # Contraction     
        elif r_fitness >= population[-2].fitness:
            # Outside Contraction
            if r_fitness < population[-1].fitness:
                oc_contract = centroid + d_oc * ( centroid - population[-1].chromosome)
                _, oc_fitness, oc_path = run_md(Model(oc_contract, 0))
                if oc_fitness <= r_fitness:
                    population.pop(-1)
                    population.append(Model(oc_contract, oc_fitness))
                    simplex_history.append((oc_fitness, oc_path))
                    with open("pop_history.log", "a") as fout:
                        print("Outside contraction  ("+str(oc_contract)+"), fitness: "+str(oc_fitness), file=fout)

            # Inside Contraction
            elif r_fitness >= population[-1].fitness:
                ic_contract = centroid + d_ic * ( centroid - population[-1].chromosome)
                _, ic_fitness, ic_path = run_md(Model(ic_contract, 0))
                if ic_fitness <= r_fitness:
                    population.pop(-1)
                    population.append(Model(ic_contract, ic_fitness))
                    simplex_history.append((ic_fitness, ic_path))
                    with open("pop_history.log", "a") as fout:
                        print("Inside contraction ("+str(ic_contract)+"), fitness: "+str(ic_fitness), file=fout)

                else: # Shrink
                    with open("pop_history.log", "a") as fout:
                        print("Simplex shrinking, contraction ("+str(ic_contract)+"), fitness: "+str(ic_fitness), file=fout)
                    for i, j in enumerate(population[1:]):
                        new_chromosome = population[0].chromosome + d_s * (np.subtract(j.chromosome, population[0].chromosome))
                        _, new_fitness, new_path = run_md(Model(new_chromosome, 0))
                        population.pop(1)
                        population.append(Model(new_chromosome, new_fitness))
                        simplex_history.append((new_fitness, new_path))

    with open("pop_history.log", "a") as fout:
        print("\nExplored models: ", file=fout)
        print("-------------------------------------------------------------------------------", file=fout)
        simplex_history.sort()
        for model in simplex_history:
            print(model[0],">>>",model[1], file=fout)
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", file=fout)

# gen_alg(population)

# # if argv -l locate and read pop_history.log
population = Model.from_log(input)

# # if various pop_history.log in dir
#population = [model for root, dirs, files in os.walk(input) for file in files if file.endswith(".log") and file.startswith("pop_history") for model in Model.from_log(os.path.join(root, file))]

# # else look for xvgs
# population = [Model.from_xvg(os.path.join(root, file)) for root, dirs, files in os.walk(input) for file in files if file.endswith(".xvg") and file.startswith("condensation_")]

simplex_nelder_mead(population, n_vars=0, n_iter=20)

exit(1)

import scipy as sp

sp.optimize.minimize()
