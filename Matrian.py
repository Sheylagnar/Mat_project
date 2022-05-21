import csv
import os
from sys import platform as sys_pf
if sys_pf == "darwin":
    import matplotlib
    matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Bio import SeqIO
from math import log, sqrt
from collections import defaultdict
import logging

class matrian:  # lower matrix
    """calculate and print main genetic distances results.
    Parameters
    ----------
    path: str
        The path to folder with the input file.
    fasta: str
        The name of fasta file.
    gen: int
        genus position on sequence names splited by "_".
    sp: int
        genus position on sequence names splited by "_".
    distance: str
        Substitution model, k for K2p or p for p-distance (default=k).
    """

    def __init__(self, path, fasta, gen, sp, distance, out_name):
        self.path = path
        self.fasta = open(path + fasta, "r")
        self.out_name = out_name
        self.fasta_seqIO = SeqIO.parse(self.fasta, "fasta")  # sequence object
        self.fasta_names = [
            seq.id for seq in self.fasta_seqIO]
        self.lst_tot = self.lista_total()
        self.ncol = self.lst_tot[0][1]
        self.fasta_dict = SeqIO.to_dict(
            SeqIO.parse(fasta, "fasta"))
        self.Lname = self.name_sp(gen, sp)
        self.Lsp = self.listar_sp(sp)
        self.Lgen = self.listar_sp(gen)
        self.data = self.matrix(distance)
        self.data2 = self.matrix_new(distance)

    def fasta_seq(self):  # genera lista de tuplas de nombre y sequencia
        self.fasta.seek(0)
        aln = self.fasta
        name = []
        seq = []
        D = []
        for line in aln:
            line = line.strip()
            if line.startswith(">"):
                name.append(line[1:])
                D += seq
                seq = ['']
            else:
                if not line:
                    break
                if line[0] != ">":
                    seq[0] += (line.rstrip())
        D += seq  # apendea la ultima secuencia
        A = list(zip(name, D))
        return A

    def K2Pdistance(self, seq1,
                    seq2):  # adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
        pairs = []
        for x in zip(seq1, seq2):
            if "-" not in x and "N" not in x:
                pairs.append(x)
        ts_count = 0
        tv_count = 0
        length = len(pairs)
        transitions = ["AG", "GA", "CT", "TC"]
        transversions = ["AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"]
        for (x, y) in pairs:
            if x + y in transitions:
                ts_count += 1
            elif x + y in transversions:
                tv_count += 1
        p = float(ts_count) / length
        q = float(tv_count) / length
        try:
            d = -0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q)) * 100
        except ValueError:
            logging.info("Tried to take log of a negative number")
            return None
        return d

    def pdistance(self, seq1,
                  seq2):  # adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
        p = 0
        pairs = []
        for x in zip(seq1, seq2):
            if "-" not in x and "N" not in x:
                pairs.append(x)
        for (x, y) in pairs:
            if x != y:
                p += 1
        length = len(pairs)
        return (float(p) / length) * 100

    def matrix(self, distance):  # genera la matriz de distancia
        if distance == 'k':
            print('Using k2p distance\n')
        else:
            print('Using p-distance\n')
        m = [[None] * self.ncol for i in range(self.ncol)]
        A = self.fasta_seq()
        for i in range(self.ncol):
            for j in range(self.ncol):
                #                if i==j or i>j:    ###Verificar si funciona!! y mejorar el tiempo
                #                    continue
                if distance == 'k':
                    num = self.K2Pdistance(A[i][1], A[j][1])
                    if num == -0.0:
                        num = 0.0
                    m[i][j] = (round(num, 5))  # code by neguinha!
                if distance == 'p':
                    num = self.pdistance(A[i][1], A[j][1])
                    if num == -0.0:
                        num = 0.0
                    m[i][j] = (round(num, 5))  # code by neguinha!
        return m

    def matrix_new(self, distance):
        if distance == "k":
            logging.info("Using k2p distance\n")
        else:
            logging.info("Using p-distance\n")
        dmatrix = defaultdict(dict)
        for i in self.fasta_names:
            for j in self.fasta_names:
                if distance == "k":
                    num = self.k2Pdistance(
                        self.fasta_dict[i].seq, self.fasta_dict[j].seq
                    )
                    if num == -0.0:
                        num = 0.0
                    dmatrix[i][j] = round(num, 5)  # code by neguinha!
                if distance == "p":
                    num = self.pdistance(
                        self.fasta_dict[i].seq, self.fasta_dict[j].seq
                    )
                    if num == -0.0:
                        num = 0.0
                    dmatrix[i][j] = round(num, 5)  # code by neguinha!
        return dmatrix

    def padron(self, i):  # genera lista con padron '_'
        Padron = [0]
        Padro_pos = 0
        for letra in i:
            Padro_pos += 1
            if letra == '_':
                Padron.append(Padro_pos)
        Padron.append(self.ncol)
        return Padron

    def listar_sp(self, a):  # genera lista de especies o generos
        self.fasta.seek(0)
        L = []
        sp = self.fasta_names
        #        print sp
        for i in sp:
            Padron = self.padron(i)
            ini = (Padron[a - 1] + 1)
            fin = (Padron[a] - 1)
            L.append(i[ini - 1:fin])
        i1 = ''
        pos = 0
        L1 = []
        for i in (L):
            pos += 1
            if i != i1:
                L1.append(pos)
                i1 = L[pos - 1]
        L2 = []
        for i in range(len(L1) - 1):
            L3 = [L1[i], L1[i + 1] - 1]
            L2.append(L3)
        L2.append([L1[len(L1) - 1], len(L)])
        return L2

    def lista_total(self):
        self.fasta.seek(0)
        T = self.fasta_names
        L = [[1, len(T)]]
        return L

    def name_sp(self, a, b):  # genera lista de nombres de especies
        self.fasta.seek(0)
        L = []
        sp = self.fasta_names
        for i in sp:
            Padron = self.padron(i)
            ini = (Padron[a - 1] + 1)
            fin2 = (Padron[b] - 1)
            L.append(i[ini - 1:fin2])
        Lname = []
        [Lname.append(key) for key in L if key not in Lname]
        return Lname

    def conc_intra(self, ini, fin):  # concatena intra dentro de una sp
        if fin > ini:
            x = ini - 1
            y = ini - 1
            z = ini - 1
            conc = []
            while x < fin and z < fin:
                while y <= z and y != x:
                    conc.append(self.data[x][y])
                    y += 1
                z += 1
                y = ini - 1
                x += 1
            return conc
        else:
            return None

    def conc_inter(self, ini, fin, ini2, fin2):  # conc inter de um sp por genero
        if ini == ini2 and fin == fin2:
            return None
        x = ini - 1
        y = ini2 - 1
        conc = []
        while x < fin:
            while y < ini - 1:
                conc.append(self.data[x][y])
                y += 1
            y = ini2 - 1
            x += 1
        x = fin
        y = ini - 1
        while x < fin2:
            while y < fin:
                conc.append(self.data[x][y])
                y += 1
            y = ini - 1
            x += 1
        return conc

    def intra(self, ini, fin):  # da la media intra para un dado grupo
        L = self.conc_intra(ini, fin)
        if L == None:
            return None
        v = 0
        for i in L:
            v += i
        m = (float(v)) / (len(L))
        return m

    def cont_intra_all(self):  # cont freq intra de todos los grupos
        cont = []
        for i in self.Lsp:
            if i[1] > i[0]:
                cont += (self.conc_intra(i[0], i[1]))
        cont.sort()
        # cont2={x:cont.count(x) for x in set(cont)}
        # return cont2
        return cont

    def intra_all(self):  # da media intra normalizada de varios grupos #used
        cont_intra = self.cont_intra_all()
        M = round((sum(cont_intra) / len(cont_intra)), 5)
        return M

    def inter(self, ini, fin, ini2, fin2):  # da media interespecifica por genero
        if ini == ini2 and fin == fin2:
            return None
        L = self.conc_inter(ini, fin, ini2, fin2)
        v = 0
        for i in L:
            v += i
        m = (float(v)) / (len(L))
        return m

    def cont_inter_all(self):  # cont freq inter de todos los grupos dentro de su genero
        cont_inter = []
        for g in self.Lgen:  # Lgenlst_tot
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    if (self.conc_inter(e[0], e[1], g[0], g[1])) != None:
                        cont_inter += (self.conc_inter(e[0], e[1], g[0], g[1]))
        cont_inter.sort()
        # cont2={x:(cont_inter.count(x)/2) for x in set(cont_inter)}
        # return cont2
        return cont_inter

    def inter_all(self):  # da media intergrupo normalizada #used
        cont_inter = self.cont_inter_all()
        M = round((sum(cont_inter) / len(cont_inter)), 5)
        return M

    def interALL_all(self):  # da media inter general normalizada #used
        m = 0
        n = 0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    if self.inter(e[0], e[1], g[0], g[1]) != None:
                        m += self.inter(e[0], e[1], g[0], g[1])
                        n += 1
        M = m / n
        return M

    def inter_intra_all_name(self):  # da media inter e intra de todos los grupos y nombre #used
        II = []
        x = 0
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    name = self.Lname[x]
                    x += 1
                    inter = self.inter(e[0], e[1], g[0], g[1])
                    intra = self.intra(e[0], e[1])
                    intera = [name, inter, intra]
                    II.append(intera)
        return II

    def interALL_intra_all_name(self):  # da media inter general e intra de todos los grupos y nombre #used
        II = []
        x = 0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    name = self.Lname[x]
                    x += 1
                    inter = self.inter(e[0], e[1], g[0], g[1])
                    intra = self.intra(e[0], e[1])
                    intera = [name, inter, intra]
                    II.append(intera)
        return II

    def max_intra(self, ini, fin):  # da maxima dentro de un grupo
        L = self.conc_intra(ini, fin)
        if L == None:
            return None
        v = max(L)
        return v

    def min_intra(self, ini, fin):  # da minima dentro de un grupo
        L = self.conc_intra(ini, fin)
        if L == None:
            return None
        v = min(L)
        return v

    def min_inter(self, ini, fin, ini2, fin2):  # da minima inter de um grupo
        L = self.conc_inter(ini, fin, ini2, fin2)
        if L == None:
            return None
        v = min(L)
        return v

    def max_inter(self, ini, fin, ini2, fin2):  # da maxima inter de um grupo
        L = self.conc_inter(ini, fin, ini2, fin2)
        if L == None:
            return None
        v = max(L)
        return v

    def min_inter_pos(self, ini, fin, ini2, fin2):  # posicion en lista conc_inter del minima inter de um grupo #used
        L = self.conc_inter(ini, fin, ini2, fin2)
        if L == None:
            return None
        vpos = L.index(min(L))
        return vpos

    def min_all2(self):  # da min intergrupos #used
        MM = []
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    mi = self.min_inter(e[0], e[1], g[0], g[1])
                    if mi != None:
                        MM.append(mi)
        return MM

    def max_inter_all(self):  # da max intergrupos #used
        MM = []
        for g in self.Lgen:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    mi = self.max_inter(e[0], e[1], g[0], g[1])
                    MM.append(mi)
        return MM

    def max_inter_ALL(self):  # da max inter sin importar los grupos #used
        MM = []
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    mi = self.max_inter(e[0], e[1], g[0], g[1])
                    MM.append(mi)
        return MM

    def min_ALL(self):  # da min inter general sin importar los grupos #used
        MM = []
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    mi = self.min_inter(e[0], e[1], g[0], g[1])
                    MM.append(mi)
        return MM

    def max_all(self):  # da max intra de todos los grupos #used
        MM = []
        for i in self.Lsp:
            ma = self.max_intra(i[0], i[1])
            MM.append(ma)
        return MM

    def min_intra_all(self):  # da min intra de todos los grupos #used
        MM = []
        for i in self.Lsp:
            ma = self.min_intra(i[0], i[1])
            if ma != None:
                MM.append(ma)
        return MM

    def min_max_ALL_name(self):  # da min general y max intra con nombre de sp       #used
        MM = []
        x = 0
        for g in self.lst_tot:
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    no = self.Lname[x]
                    x += 1
                    mi = self.min_inter(e[0], e[1], g[0], g[1])
                    ma = self.max_intra(e[0], e[1])
                    mima = [no, mi, ma]
                    MM.append(mima)
        return MM

    def min_ALL_pos(self):  # da min general y max con nombre de sp      #used
        min_inter_dict = {}
        MM_name_indv_NN = {}
        for sp in self.Lname:
            ind_sp = [
                ind for ind in self.fasta_names if ind.startswith(sp + '_')
            ]  # list of individuals of a given sp
            L = []
            K = []
            for ind in ind_sp:
                ind_dic = self.data2[
                    ind
                ]  # dictionary of individuals genetic distances
                L = L + [
                    values
                    for key, values in ind_dic.items()
                    if not key.startswith(sp + '_')
                ]
                K = K + [
                    key for key in ind_dic.keys() if not key.startswith(sp + '_')
                ]
            v = min(L)
            kpos = L.index(v)
            min_inter_dict[sp] = v
            MM_name_indv_NN[sp] = K[kpos]
        return min_inter_dict, MM_name_indv_NN

    def min_media_max_intra(self):  # used
        MX = max(self.max_all())
        MI = min(self.min_intra_all())
        med = self.intra_all()
        return MI, med, MX

    def min_media_max_inter(self):  # used
        MX = max(self.max_inter_all())
        MI = min(self.min_all2())
        med = self.inter_all()
        return MI, med, MX

    def min_media_max_tot(self):  # used
        MX = max(self.max_inter_ALL())
        MI = min(self.min_ALL())
        med = self.interALL_all()
        return MI, med, MX

    def inter_geo(self):  # da valores unico inter de todos los grupos dentro de su grupo mayor
        cont_inter = []
        x = 0
        for g in self.Lgen:  # Lgenlst_tot
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    if (self.conc_inter(e[0], e[1], g[0], g[1])) != None:
                        inter = (self.conc_inter(e[0], e[1], g[0], g[1]))
                        inter.sort()
                        inter2 = []
                        for e in inter:
                            if e not in inter2:
                                inter2.append(e)
                        no = self.Lname[x]
                        x += 1
                        inno = [no, inter2]
                        cont_inter.append(inno)
        return cont_inter

    def cont_inter_geo(self):  # cont freq unico inter de todos los grupos dentro de su grupo mayor
        cont_inter = []
        for g in self.Lgen:  # Lgenlst_tot
            for e in self.Lsp:
                if e[1] <= g[1] and e[0] >= g[0]:
                    if (self.conc_inter(e[0], e[1], g[0], g[1])) != None:
                        inter = (self.conc_inter(e[0], e[1], g[0], g[1]))
                        inter.sort()
                        inter2 = []
                        for e in inter:
                            if e not in inter2:
                                inter2.append(e)
                        cont_inter.append(inter2)
        y = 0
        cont2 = []
        cont3 = []
        for e in range(len(self.Lsp)):
            if y == e:
                cont2.append(cont_inter[y])
                y += 2
        for e in cont2:
            for i in e:
                cont3.append(i)
        cont3.sort()
        cont4 = {x: (cont3.count(x)) for x in set(cont3)}
        return cont4

    def inter_geo_max(self):  # da valores unico inter de todos los grupos dentro de su grupo mayor
        max = self.max_inter_all()
        print(max)

    def plot_freq(self):  # Barcoding gap graph
        ter = self.cont_inter_all()
        tra = self.cont_intra_all()
        newBins_tra = len(set(tra)) // 3
        if newBins_tra == 0:
            newBins_tra = 1
        newBins_ter = len(set(ter)) // 3
        f, ax = plt.subplots(3, 1, sharex='col', sharey='all')
        sns.histplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance', ax=ax[1])
        sns.histplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance', ax=ax[0])
        sns.histplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance', ax=ax[2])
        sns.histplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance', ax=ax[2])
        ax[0].set_title('DNA Barcoding gap')
        ax[0].set_ylabel('# of taxon pairs')
        ax[1].set_ylabel('# of taxon pairs')
        ax[2].set_ylabel('# of taxon pairs')
        ax[0].legend()
        ax[1].legend()
        ax[2].legend()
        ax[2].set_xlabel('Genetic Distance')
        f.savefig(os.path.join(self.path, self.out_name+'_barcoding_gap.pdf'))
        #plt.show(f)
        #plt.close(f)
        plt.clf()

    def med_ind_sp(self):  # media de individuso por sp
        m = 0
        t = 0
        for i in self.Lsp:
            m += i[1] - i[0] + 1
            t += 1
        return float(m / t)

    def plot_max_min(self, df):  # max vs min graph
        sns.lmplot(x='intra2', y='inter',
                   data=df,
                   fit_reg=False)
        plt.title('Maximum intraspecific vs Minimum to NN')
        plt.xlabel('Maximum intraspecific')
        plt.ylabel('Minimum to NN')
        z = [df["inter"].max(), df["intra2"].max()]
        plt.axis([0, max(z) + 1, 0, max(z) + 1])
        lims = [0, max(z) + 1]
        plt.plot(lims, lims, ':k')
        plt.savefig(os.path.join(self.path, self.out_name + '_min_max.pdf'))
        # plt.show()
        plt.clf()

    def analyze(self):  # hace analisis en bloque
        print('Summary table (Name, mean intra, max intra, NN, distance to NN) in percentage')
        a, b, c, d, e = [], [], [], [], []
        II = self.inter_intra_all_name()
        mima = self.min_max_ALL_name()
        name_NNindv = self.min_ALL_pos()
        for i in range(len(II)):
            a.append(II[i][0])
            b.append(II[i][2])
            c.append(mima[i][2])
            e.append(mima[i][1])
        for i in self.Lname:
            d.append(name_NNindv[1].get(i))
        data = {"Name": a, "Mean": b, "Max": c, "NN": d, "DtoNN": e}
        summ = pd.DataFrame(data)
        print(summ.to_string(index=False))
        print('')
        print('Min interspecific and max intraspecific by group')
        mima = self.min_max_ALL_name()
        inter = []
        intra2 = []
        intra = []
        name = []
        for i in mima:
            inter.append(i[1])
            intra.append(i[2])
            if i[2] == None:
                intra2.append(0.0)
            else:
                intra2.append(i[2])
            name.append(i[0])
        datas2 = {"name": name, "inter": inter, "intra": intra, "intra2": intra2}
        df = pd.DataFrame(datas2)
        df1 = df[['name', 'inter', 'intra']].copy()
        print(df1.to_string(index=False))
        print("")
        print('Min, max and mean at all')
        tra = []
        ter = []
        title = ['minimum', 'mean', 'maximum']
        L = self.min_media_max_intra()
        for i in L:
            tra.append(i)
        M = self.min_media_max_inter()
        for i in M:
            ter.append(i)
        datas3 = {"intra": tra, "inter": ter}
        tab = pd.DataFrame(datas3, index=title)
        print(tab.transpose())
        ####Plot max vc min graph####
        self.plot_max_min(df)
        ####Plot frequencies graph####
        self.plot_freq()

    def geoanalyze(self):  # hace analisis geobarcode
        print('valores inter')
        ter = self.inter_geo()
        for i in ter:
            print(i)

    #        ter2=self.cont_inter_geo()
    #        for k,v in ter2.items():
    #            print k,v

    def geoanalyze_intra(self):  # hace analisis geobarcode
        print('valores intra BIN')
        ter2 = self.cont_intra_all()
        cont = ter2
        cont.sort()
        cont2 = {x: (cont.count(x)) for x in set(cont)}
        sumv = 0
        for k, v in cont2.items():
            sumv += v
        x = []
        y = []
        z = []
        for k, v in cont2.items():
            x.append(k)
            y.append(v * 100 / float(sumv))
            z.append(v)
        data = {"intra": x, "relative_freq": y, 'absolute_fre': z}
        df = pd.DataFrame(data)
        df.to_csv(self.out_name+"_intra.csv")
        print("done")
        return df


    def geoanalyze_max_intra(self):  # hace analisis geobarcode
        print('valores max intra BIN')
        mima = self.min_max_ALL_name()
        x = []
        y = []
        y2 = []
        z = []
        for i in mima:
            x.append(i[1])
            y2.append(i[2])
            if i[2] == None:
                y.append(0.0)
            else:
                y.append(i[2])
            z.append(i[0])
        cont = y2
        cont.sort()
        cont2 = {x: (cont.count(x)) for x in set(cont)}
        sumv = 0
        for k, v in cont2.items():
            sumv += v
        x = []
        y = []
        z = []
        for k, v in cont2.items():
            x.append(k)
            y.append(v * 100 / float(sumv))
            z.append(v)
        #        df=df.astype(float)
        data = {"intra": x, "relative_freq": y, 'absolute_fre': z}
        df = pd.DataFrame(data)
        #        print df1.to_string(index=False)
        df.to_csv(self.out_name + "_max_intra.csv")
        return df


#    def geoanalyze_inter(self):
#        self.inter_geo_max()
#        m=self.max_all()
#        print m
#        print self.Lsp

def main(path, fasta, gen, sp, distance, n=False, geo_m=False, out_name=False):
    tmp = matrian(path, fasta, gen, sp, distance, out_name)
#    tmp.analyze()
#     tmp.geoanalyze_inter()
#     tmp.geoanalyze()
    if geo_m == 'i':
        df = tmp.geoanalyze_intra()
        return df
    if geo_m == 'j':
        df = tmp.geoanalyze_max_intra()
        return df


#    if n==True:
#        List=tmp.name_sp(gen,sp)
#        with open(out_name,'w+') as nominal_file:
#            for i in List:
#                nominal_file.write(i+'\n')
if __name__ == "__main__":
    main()