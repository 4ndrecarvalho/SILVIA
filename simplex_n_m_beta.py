#! /usr/bin/env python3
import numpy as np
import sys
import os

input=sys.argv[1]
# directory = "/home/acarvalho/Desktop/silvia/md/stella8VS_new_topology_Dall/"

#Valores ~ experimentais de Devreux et al.
max_c = .9
d_q_max_q0 = 1
d_q_max_q1 = .58
d_q_max_q2 = .58
d_q_max_q3 = .60
d_q_max_q4 = .47
d_q_cross_q0_q1 = .45
d_q_cross_q1_q2 = .45
d_q_cross_q2_q3 = .44
d_q_cross_q3_q4 = .50
d_q_cross_q0_q2 = .20
d_q_cross_q1_q3 = .22
d_q_cross_q2_q4 = .22
d_q_cross_q0_q3 = .06
d_q_cross_q1_q4 = .07
d_q_cross_q0_q4 = .01
d_c_max_q1 = .25
d_c_max_q2 = .52
d_c_max_q3 = .77
d_c_max_q4 = 0.90
d_c_cross_q0_q1 = .17
d_c_cross_q1_q2 = .41
d_c_cross_q2_q3 = .64
d_c_cross_q3_q4 = .87
d_c_cross_q0_q2 = .26
d_c_cross_q1_q3 = .32
d_c_cross_q2_q4 = .75
d_c_cross_q0_q3 = .36
d_c_cross_q1_q4 = .62
d_c_cross_q0_q4 = .46

def crossing_points_list(qi, qj, c, cutoff, name):
    """
    Avalia os pontos de cruzamento entre listas
    cutoffs foram aplicados para remover cruzamentos
    inuteis, como os das formaçoes iniciais de q3 e q4
    ou os cruzamentos finais entre especies q0 e q1.
    """
    q_value = c_index = c_cross_qi_qj = 0
    qi=np.array(qi)
    qj=np.array(qj)
    qi[qi<cutoff] = 0.001
    qj[qj<cutoff] = 0.0001
    diff=(qj-qi)
    indexes = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]+1 #indexes where sign changes
    if len(indexes) < 1:
        q_value = 0 
        c_cross_qi_qj = 0
        # print("No intersection in",name)
    else:
        c_index = int(np.round(np.mean(indexes)))
        c_cross_qi_qj = c[c_index]
        for i in indexes:
            q_value += (qi[i] + qj[i] + qi[i-1] + qj[i-1])/4
        q_value = q_value/len(indexes)
    # print("{:s} >> q_value: {:.2f}%, c_cross_qi_qj: {:.2f}%, indexes: {:d}, c_index: {:d}".format(name,q_value*100,c_cross_qi_qj*100,len(indexes),c_index))
    return q_value, c_cross_qi_qj

def xvg_reader(fname):
    """
    Leitor de xvg de condensação obtem informação
    sobre os máximos de condensação de cada espécio de Qn
    e os pontos de cruzamento entre espécies.
    """
    # print(fname)
    t  = []; q0 = []; q1 = []; q2 = []; q3 = []; q4 = []; qx = []; c  = []
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
        q_max_q0 = np.max(q0)
        q_max_q1 = np.max(q1)
        q_max_q2 = np.max(q2)
        q_max_q3 = np.max(q3)
        q_max_q4 = np.max(q4)
        q_max_qx = np.max(qx)
        c_max_q0 = c[q0.index(q_max_q0)]
        c_max_q1 = c[q1.index(q_max_q1)]
        c_max_q2 = c[q2.index(q_max_q2)]
        c_max_q3 = c[q3.index(q_max_q3)]
        c_max_q4 = c[q4.index(q_max_q4)]
        c_max_qx = c[qx.index(q_max_qx)]
        q_cross_q0_q1, c_cross_q0_q1 = crossing_points_list(q0, q1, c, .00, "q0_q1")
        q_cross_q1_q2, c_cross_q1_q2 = crossing_points_list(q1, q2, c, .00, "q1_q2")
        q_cross_q2_q3, c_cross_q2_q3 = crossing_points_list(q2, q3, c, .00, "q2_q3")
        q_cross_q3_q4, c_cross_q3_q4 = crossing_points_list(q3, q4, c, .04, "q3_q4") #filtrar formações iniciais
        q_cross_q0_q2, c_cross_q0_q2 = crossing_points_list(q0, q2, c, .00, "q0_q2")
        q_cross_q1_q3, c_cross_q1_q3 = crossing_points_list(q1, q3, c, .00, "q1_q3")
        q_cross_q2_q4, c_cross_q2_q4 = crossing_points_list(q2, q4, c, .00, "q2_q4")
        q_cross_q0_q3, c_cross_q0_q3 = crossing_points_list(q0, q3, c, .00, "q0_q3")
        q_cross_q1_q4, c_cross_q1_q4 = crossing_points_list(q1, q4, c, .00, "q1_q4")
        q_cross_q0_q4, c_cross_q0_q4 = crossing_points_list(q0, q4, c, .00, "q0_q4")
    return max_c, q_max_q0, q_max_q1, q_max_q2, q_max_q3, q_max_q4, c_max_q0, c_max_q1, c_max_q2, c_max_q3, c_max_q4, c_max_qx, q_cross_q0_q1, c_cross_q0_q1, q_cross_q1_q2, c_cross_q1_q2, q_cross_q2_q3, c_cross_q2_q3, q_cross_q3_q4, c_cross_q3_q4, q_cross_q0_q2, c_cross_q0_q2, q_cross_q1_q3, c_cross_q1_q3, q_cross_q2_q4, c_cross_q2_q4, q_cross_q0_q3, c_cross_q0_q3, q_cross_q1_q4, c_cross_q1_q4, q_cross_q0_q4, c_cross_q0_q4, q_max_qx

def fitness(fname, a=.25, b=.25, c=.25, d=.25, e=5.):
    """ 
    Função de aptidão é baseada nos valores de
    condensação de cada espécie em função do grau
    de condensação(c) e da fracção molar (q).
    Máximos de cada espécie, pontos de cruzamento
    entre espécies.
    """
    max_c, q_max_q0, q_max_q1, q_max_q2, q_max_q3, q_max_q4, c_max_q0, c_max_q1, c_max_q2, c_max_q3, c_max_q4, c_max_qx, q_cross_q0_q1, c_cross_q0_q1, q_cross_q1_q2, c_cross_q1_q2, q_cross_q2_q3, c_cross_q2_q3, q_cross_q3_q4, c_cross_q3_q4, q_cross_q0_q2, c_cross_q0_q2, q_cross_q1_q3, c_cross_q1_q3, q_cross_q2_q4, c_cross_q2_q4, q_cross_q0_q3, c_cross_q0_q3, q_cross_q1_q4, c_cross_q1_q4, q_cross_q0_q4, c_cross_q0_q4, q_max_qx= xvg_reader(fname)

    q_max_qs = .9*((np.absolute(q_max_q1-d_q_max_q1)/d_q_max_q1) + (np.absolute(q_max_q2-d_q_max_q2)/d_q_max_q2) + (np.absolute(q_max_q3-d_q_max_q3)/d_q_max_q3)) + 0.1*(np.absolute(q_max_q4-d_q_max_q4)/d_q_max_q4)
    c_max_qs = .9*((np.absolute(c_max_q1-d_c_max_q1)/d_c_max_q1) + (np.absolute(c_max_q2-d_c_max_q2)/d_c_max_q2) + (np.absolute(c_max_q3-d_c_max_q3)/d_c_max_q3)) + 0.1*(np.absolute(c_max_q4-d_c_max_q4)/d_c_max_q4)
    q_cross_points_top = (1/4)*((np.absolute(q_cross_q1_q2-d_q_cross_q1_q2)/d_q_cross_q1_q2) + (np.absolute(q_cross_q2_q3-d_q_cross_q2_q3)/d_q_cross_q2_q3)) + (3/4)*((np.absolute(q_cross_q3_q4-d_q_cross_q3_q4)/d_q_cross_q3_q4))
    q_cross_points_middle = (1/3)*((np.absolute(q_cross_q0_q2-d_q_cross_q0_q2)/d_q_cross_q0_q2) + (np.absolute(q_cross_q1_q3-d_q_cross_q1_q3)/d_q_cross_q1_q3) + (np.absolute(q_cross_q2_q4-d_q_cross_q2_q4)/d_q_cross_q2_q4))
    q_cross_points_bottom = .5*((np.absolute(q_cross_q0_q3-d_q_cross_q0_q3)/d_q_cross_q0_q3) + (np.absolute(q_cross_q1_q4-d_q_cross_q1_q4)/d_q_cross_q1_q4))
    c_cross_points_top = (1/4)*((np.absolute(c_cross_q1_q2-d_c_cross_q1_q2)/d_c_cross_q1_q2) + (np.absolute(c_cross_q2_q3-d_c_cross_q2_q3)/d_c_cross_q2_q3)) + (3/4)*((np.absolute(c_cross_q3_q4-d_c_cross_q3_q4)/d_c_cross_q3_q4))
    c_cross_points_middle = (1/3)*((np.absolute(c_cross_q0_q2-d_c_cross_q0_q2)/d_c_cross_q0_q2) + (np.absolute(c_cross_q1_q3-d_c_cross_q1_q3)/d_c_cross_q1_q3) + (np.absolute(c_cross_q2_q4-d_c_cross_q2_q4)/d_c_cross_q2_q4))
    c_cross_points_bottom = .5*((np.absolute(c_cross_q0_q3-d_c_cross_q0_q3)/d_c_cross_q0_q3) + (np.absolute(c_cross_q1_q4-d_c_cross_q1_q4)/d_c_cross_q1_q4))

    max_qs = .5*(c_max_qs + q_max_qs)
    cross_points_top = .5*(q_cross_points_top + c_cross_points_top)
    cross_points_middle = .5*(q_cross_points_middle + c_cross_points_middle)
    cross_points_bottom = .5*(q_cross_points_bottom + c_cross_points_bottom)
    
    fitness = 1 - (a*max_qs + b*cross_points_top + c*cross_points_middle + d*cross_points_bottom + e*q_max_qx) 
    # print(fname)
    # print("max_c:",max_c,"\nq_max_q0:",q_max_q0,"\nq_max_q1:",q_max_q1,"\nq_max_q2:",q_max_q2,"\nq_max_q3:",q_max_q3,"\nq_max_q4:", q_max_q4,"\nq_max_qx:", q_max_qx,"\n")
    # print("c_max_q0:", c_max_q0,"\nc_max_q1:", c_max_q1,"\nc_max_q2:", c_max_q2,"\nc_max_q3:", c_max_q3,"\nc_max_q4:", c_max_q4,"\nc_max_qx:", c_max_qx,"\n")
    # print("q_cross_q0_q1:",q_cross_q0_q1,"\nc_cross_q0_q1:",c_cross_q0_q1,"\nq_cross_q1_q2:",q_cross_q1_q2,"\nc_cross_q1_q2:",c_cross_q1_q2,"\nq_cross_q2_q3:",q_cross_q2_q3,"\nc_cross_q2_q3:",c_cross_q2_q3,"\nq_cross_q3_q4:",q_cross_q3_q4,"\nc_cross_q3_q4:",c_cross_q3_q4,"\n")
    # print("q_cross_q0_q2:",q_cross_q0_q2,"\nc_cross_q0_q2:",c_cross_q0_q2,"\nq_cross_q1_q3:",q_cross_q1_q3,"\nc_cross_q1_q3:",c_cross_q1_q3,"\nq_cross_q2_q4:",q_cross_q2_q4,"\nc_cross_q2_q4:",c_cross_q2_q4,"\n")
    # print("q_cross_q0_q3:",q_cross_q0_q3,"\nc_cross_q0_q3:",c_cross_q0_q3,"\nq_cross_q1_q4:",q_cross_q1_q4,"\nc_cross_q1_q4:",c_cross_q1_q4,"\n")
    # print("q_cross_q0_q4:",q_cross_q0_q4,"\nc_cross_q0_q4:",c_cross_q0_q4)
    # print("max_qs:",max_qs,"\ncross_points_top:",cross_points_top,"\ncross_points_middle:",cross_points_middle,"\ncross_points_bottom:",cross_points_bottom)
    # print("Fitness:",f*100,"\n")
    # print(fname,f*100)
    return fitness

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

def crossover(parent1, parent2):
    xoverpt = np.random.randint(0, 4)
    p1a,p1b = parent1.chromosome[:xoverpt], parent1.chromosome[xoverpt:]
    p2a,p2b = parent2.chromosome[:xoverpt], parent2.chromosome[xoverpt:]
    child1 = Model(p1a+p2b, np.inf)
    child2 = Model(p2a+p1b, np.inf)
    return child1, child2

def mutation(child, pmut=.5):
    for i,(_,(xmin,xmax)) in enumerate(zip(child.chromosome, child.constraints)):
        if np.random.uniform(0,1) < pmut:
            new = child.chromosome[i] + np.random.normal(0, spread)
            while new < xmin or new > xmax:
                new = child.chromosome[i] + np.random.normal(0, spread)
            child.chromosome[i] = new
    return child

def run_md(child, directory="/home/acarvalho/Desktop/silvia/simplex/"):
    path = directory+"dVS_"+str(child.chromosome[1])+"/"+str(child.chromosome[0])+"_"+str(child.chromosome[2])+"/"
    # os.makedirs(path, exist_ok=True)
    os.system("bash "+directory+"run "+str(child.chromosome[0])+" "+str(child.chromosome[2])+" "+str(child.chromosome[3])+" "+str(child.chromosome[1]))
    os.system("gmx -nocopyright -nobackup trjconv -f "+path+str(child.chromosome[3])+".xtc -s "+path+str(child.chromosome[3])+".tpr -o "+path+"d.xtc -pbc mol -skip 2 -n "+directory+"index.ndx")
    os.system("python3 "+directory+"condensation.py "+path+"d.xtc 5.8 no "+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs") #Check if Qsi=5.8, SQsi=5.2, XQsi=3.2
    os.system("mv conde*rep* "+path)
    try:
        child.fitness = fitness(path+"condensation_"+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs.xvg")
    except:
        child.fitness = -999
    os.system("echo Fitness: "+str(child.fitness)+" >> "+path+"gene_dVS_"+str(child.chromosome[1])+"_"+str(child.chromosome[0])+"_"+str(child.chromosome[2])+"_"+str(child.chromosome[3])+".log")
    child.path = path+"condensation_"+str(child.chromosome[3])+"rep_"+str(child.chromosome[1])+"dvs.xvg"
    return child.chromosome, child.fitness, child.path

def gen_alg(population, n_generations=10, n_pairs=6, n_survivors=3):
    for generation in range(n_generations):
        population.sort(key=lambda i: i.fitness, reverse=True)
        parents1 = population[0:2*n_pairs:2]
        parents2 = population[1:2*n_pairs:2]
        children = []
        for p1, p2 in zip(parents1, parents2):
            child1,child2 = crossover(p1,p2)
            mutation(child1)
            mutation(child2)
            child1.generation=child2.generation=generation+1
            children.append(child1)
            children.append(child2)
        print([(i.chromosome, i.fitness, i.generation) for i in children])
        for child in children:
            run_md(child, directory)
        population = population[:n_survivors] + children

def simplex_nelder_mead(population, n_iter=10, n_models=0, d_r=1, d_e=2, d_oc=.5, d_ic=-.5, d_s=.5):
    population.sort(key=lambda i: i.fitness, reverse=True)
    with open("pop_history.log", "a") as fout:
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", file=fout)
        print("Initial Population of "+str(input)+": ", file=fout)
        for model in population:
            print(str(model.chromosome)+"\t>->\t"+str(model.fitness), file=fout)
    simplex_history = []
    if n_models == 0:
    # automatizar o n_models através da comparação de valores entre modelos.
        for i in population:
            for j in population:
                diff = len(np.where((np.array(i.chromosome)-np.array(j.chromosome))!=0)[0])
                if diff > n_models:
                    n_models = diff
        n_models += 1
    with open("pop_history.log", "a") as fout:
        print("\nNumber of simplex points: "+str(n_models)+"\nNumber of iterations: "+str(n_iter), file=fout)
    for iter in range(0, n_iter):
        # Order
        population.sort(key=lambda i: i.fitness, reverse=True)
        population = population[:n_models]
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
        if r_fitness <= population[0].fitness and r_fitness > population[-1].fitness:
            population.pop(-1)
            population.append(Model(reflect, r_fitness))
            with open("pop_history.log", "a") as fout:
                print("Reflected ("+str(reflect)+"), fitness: "+str(r_fitness), file=fout)

        # Expand
        elif r_fitness > population[0].fitness:
            expand = centroid + d_e * (centroid - population[-1].chromosome)
            _, e_fitness, e_path = run_md(Model(expand, 0))
            if e_fitness >= r_fitness:
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
        elif r_fitness <= population[-2].fitness:
            # Outside Contraction
            if r_fitness > population[-1].fitness:
                oc_contract = centroid + d_oc * ( centroid - population[-1].chromosome)
                _, oc_fitness, oc_path = run_md(Model(oc_contract, 0))
                if oc_fitness >= r_fitness:
                    population.pop(-1)
                    population.append(Model(oc_contract, oc_fitness))
                    simplex_history.append((oc_fitness, oc_path))
                    with open("pop_history.log", "a") as fout:
                        print("Outside contraction  ("+str(oc_contract)+"), fitness: "+str(oc_fitness), file=fout)

            # Inside Contraction
            elif r_fitness <= population[-1].fitness:
                ic_contract = centroid + d_ic * ( centroid - population[-1].chromosome)
                _, ic_fitness, ic_path = run_md(Model(ic_contract, 0))
                if ic_fitness >= r_fitness:
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

simplex_nelder_mead(population, n_models=0, n_iter=20)

exit(1)

import scipy as sp

sp.optimize.minimize()
