from Bio import SeqIO
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import subprocess
import Matrian
import numpy as np

class sorting_fasta:
    def __init__(self, fasta):
        self.fasta = open(fasta, "r")
        self.path = os.path.dirname(fasta)
        self.name = os.path.basename(fasta)
        self.outname = os.path.join(self.path, self.name + '_sorted.fasta')

    def check_fasta(self):
        print('checking if sequences are aligned')
        self.fasta.seek(0)
        handle = self.fasta
        n = len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
        print(n)
        for seq in SeqIO.parse(handle, "fasta"):
            if (len(seq.seq)) == n:
                #                print "OK"
                continue
            else:
                print(seq)
                return False

    def describe_fasta(self):
        self.fasta.seek(0)
        handle = self.fasta
        n = len(next(SeqIO.parse(handle, "fasta")))
        handle.seek(0)
        n_seq = 0
        for seq in SeqIO.parse(handle, "fasta"):
            n_seq += 1
        print('Fasta file with', n_seq, 'sequences and', n, 'base pairs')

    def sort_fasta(self):
        print('sorting fasta')
        self.fasta.seek(0)
        handle = self.fasta
        l = SeqIO.parse(handle, "fasta")
        sortedList = [f for f in sorted(l, key=lambda x: x.id)]
        file_sorted = open(self.outname, 'w+')
        for s in sortedList:
            #            print s.description
            file_sorted.write('>' + s.description + '\n')
            #            return s.description
            #            print str(s.seq)
            file_sorted.write(str(s.seq) + '\n')
        print('checking if sorting was ok')
        handle.seek(0)
        file_sorted.seek(0)
        n_aln = 0
        n_aln_sorted = 0
        for seq in SeqIO.parse(handle, "fasta"):
            n_aln += 1
        #        print n_aln
        for seq in SeqIO.parse(file_sorted, "fasta"):
            n_aln_sorted += 1
        #        print n_aln_sorted
        if n_aln == n_aln_sorted:
            print('Fasta file sorted successfuly')
            pass
        else:
            self.sort_fasta(self)
        #        handle.close()
        file_sorted.close()
        return file_sorted, self.outname

    def check_pattern(self):
        self.fasta.seek(0)
        for line in self.fasta:
            if line.startswith('>'):
                pattern = line.split('_')
                if len(pattern) < 3:
                    return False

def plot_intra_perc(data, title, axis, out):
    sns.lmplot('intra', 'freq',
               data=data,
               fit_reg=False)
    plt.title(title)
    plt.xlabel('Maximum intraspecific')
    plt.ylabel('Frecuency')
    plt.axis(axis)
    plt.savefig(out)
    plt.show()
    plt.close()

def plot_intra_perc_violin(data_load, out, X, Y, Hue):
    sns.set(style="whitegrid", palette="pastel", color_codes=True)
    sns_plot = sns.violinplot(x=X, y=Y, hue=Hue,
                              split=True, inner="quartile", cut=0,
                              palette={"Yes": "y", "No": "b"},
                              data=data_load)
    #    sns.despine(left=True)
    fig = sns_plot.get_figure()
    fig.savefig(out)
    # fig.show()
    fig.clf()

def plot_intra_boxplot(data_load, out, X, Y, Hue=False):
    sns_plot = sns.boxplot(x=X, y=Y, data=data_load)
    #    sns.despine(left=True)
    fig = sns_plot.get_figure()
    fig.savefig(out)
    # fig.show()
    fig.clf()

def plot_freq(self):  # Barcoding gap graph
    ter = self.cont_inter_all()
    tra = self.cont_intra_all()
    newBins_tra = len(set(tra)) / 3
    newBins_ter = len(set(ter)) / 3
    f, ax = plt.subplots(3, 1, sharex='col', sharey='all')
    sns.distplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance', ax=ax[1])
    sns.distplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance', ax=ax[0])
    sns.distplot(tra, bins=newBins_tra, color="r", kde=False, label='Intragroup distance', ax=ax[2])
    sns.distplot(ter, bins=newBins_ter, color="b", kde=False, label='Intergruop distance', ax=ax[2])
    ax[0].set_title('DNA Barcoding gap')
    ax[0].set_ylabel('# of taxon pairs')
    ax[1].set_ylabel('# of taxon pairs')
    ax[2].set_ylabel('# of taxon pairs')
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    ax[2].set_xlabel('Genetic Distance')
    f.savefig(self.path + 'barcoding_gap.pdf')
    #    plt.show(f)
    #    plt.close(f)
    plt.clf()

fasta='/home/shey/Escritorio/mat_project/new_db/seq_align_last.fas'
path=os.path.dirname(fasta)
fasta_sort=sorting_fasta(fasta)   #Checking fasta
fasta_sort.describe_fasta()
if fasta_sort.check_fasta() is False:
    sys.exit("Error: different length sequences, please check your alignment")
name_sorted = fasta_sort.sort_fasta()

with open(name_sorted[1]) as fasta, \
        open(os.path.join(path,'ALN_nodisp.fasta'), 'w') as newfasta_nodisp, \
        open(os.path.join(path,'ALN_paleo.fasta'), 'w')  as newfasta_paleo, \
        open(os.path.join(path,'ALN_head.fasta'), 'w') as newfasta_head, \
        open(os.path.join(path,'ALN_disp.fasta'), 'w') as newfasta_disp:
    fasta.seek(0)
    s=SeqIO.parse(fasta, "fasta")
    paleo=[]
    nodisp=[]
    head=[]
    bin_used=""
    n_pu=0
    n_bu=0
    for seq in s:
        cols=seq.id.split('_')
        if cols[1] != bin_used:
            if n_pu==1 and n_bu==1:
                nodisp.append(bin_used)
            if n_pu==1 and n_bu>1:
                paleo.append(bin_used)
            if n_pu>1 and n_bu>=1:
                head.append(bin_used)
            bin_used=cols[1]
            paleo_used=""
            n_pu=0
            basin_used=""
            n_bu=0
        else:
            if cols[2] != paleo_used:
                n_pu+=1
                paleo_used=cols[2]
            if cols[3] != basin_used:
                n_bu+=1
                basin_used=cols[3]
    if n_pu==1 and n_bu==1:
        nodisp.append(bin_used)
    if n_pu==1 and n_bu>1:
        paleo.append(bin_used)
    if n_pu>1 and n_bu>=1:
        head.append(bin_used)

    fasta.seek(0)
    s = SeqIO.parse(fasta, "fasta")
    for seq in s:
        cols = seq.id.split('_')
        for i in paleo:
            if i == cols[1]:
                newfasta_paleo.write('>' + seq.description + '\n')
                newfasta_paleo.write(str(seq.seq) + '\n')
                newfasta_disp.write('>' + seq.description + '\n')
                newfasta_disp.write(str(seq.seq) + '\n')
        for i in nodisp:
            if i == cols[1]:
                newfasta_nodisp.write('>' + seq.description + '\n')
                newfasta_nodisp.write(str(seq.seq) + '\n')
        for i in head:
            if i == cols[1]:
                newfasta_head.write('>' + seq.description + '\n')
                newfasta_head.write(str(seq.seq) + '\n')
                newfasta_disp.write('>' + seq.description + '\n')
                newfasta_disp.write(str(seq.seq) + '\n')
    fasta.seek(0)
    newfasta_nodisp.seek(0)
    newfasta_paleo.seek(0)
    newfasta_head.seek(0)
    newfasta_disp.seek(0)

fasta_head = sorting_fasta(os.path.join(path, 'ALN_head.fasta'))  # Checking fasta
if fasta_head.check_fasta() is False:
    sys.exit("Error: different length sequences, please check your alignment")
fasta_head.sort_fasta()

fasta_disp = sorting_fasta(os.path.join(path, 'ALN_disp.fasta'))  # Checking fasta
if fasta_disp.check_fasta() is False:
    sys.exit("Error: different length sequences, please check your alignment")
fasta_disp.sort_fasta()

fasta_paleo = sorting_fasta(os.path.join(path, 'ALN_paleo.fasta'))  # Checking fasta
if fasta_paleo.check_fasta() is False:
    sys.exit("Error: different length sequences, please check your alignment")
fasta_paleo.sort_fasta()

fasta_nodisp = sorting_fasta(os.path.join(path, 'ALN_nodisp.fasta'))  # Checking fasta
if fasta_nodisp.check_fasta() is False:
    sys.exit("Error: different length sequences, please check your alignment in ALN_nodisp.fasta")
fasta_nodisp.sort_fasta()

disp_df = Matrian.main(path+"/", 'ALN_disp.fasta_sorted.fasta', 1, 2, 'p', geo_m='i', out_name="disp")
nodisp_df = Matrian.main(path+"/", 'ALN_nodisp.fasta_sorted.fasta', 1, 2, 'p', geo_m='i', out_name="no_disp")
head_df = Matrian.main(path+"/", 'ALN_head.fasta_sorted.fasta', 1, 2, 'p', geo_m='i', out_name="head")
paleo_df = Matrian.main(path+"/", 'ALN_paleo.fasta_sorted.fasta', 1, 2, 'p', geo_m='i', out_name="paleo")

disp_df.to_pickle(os.path.join(path,'max_freq_disp.pkl'))
nodisp_df.to_pickle(os.path.join(path,'max_freq_nodisp.pkl'))
head_df.to_pickle(os.path.join(path,'max_freq_head.pkl'))
paleo_df.to_pickle(os.path.join(path,'max_freq_paleo.pkl'))

disp_df = pd.read_pickle(os.path.join(path,'max_freq_disp.pkl'))
nodisp_df = pd.read_pickle(os.path.join(path,'max_freq_nodisp.pkl'))
head_df = pd.read_pickle(os.path.join(path,'max_freq_head.pkl'))
paleo_df = pd.read_pickle(os.path.join(path,'max_freq_paleo.pkl'))

x=disp_df['intra']
y=disp_df['relative_freq']
z=disp_df['absolute_fre']
x2=nodisp_df['intra']
y2=nodisp_df['relative_freq']
z2=nodisp_df['absolute_fre']
x3=head_df['intra']
y3=head_df['relative_freq']
z3=head_df['absolute_fre']
x4=paleo_df['intra']
y4=paleo_df['relative_freq']
z4=paleo_df['absolute_fre']

xs = [x, x2, x3, x4]
ys = [y, y2, y3, y4]
zs = [z, z2, z3, z4]
intra = []
for i in xs:
    intra.extend(i)
freq = []
for i in ys:
    freq.extend(i)
test = []
types = []
category = []
for i in x:
    test.append('D/N')
    types.append('Yes')
    category.append('Disp')
for i in x2:
    test.append('D/N')
    types.append('No')
    category.append('No Disp')
for i in x3:
    test.append('H/P')
    types.append('Yes')
    category.append('Head')
for i in x4:
    test.append('H/P')
    types.append('No')
    category.append('Paleo')

data = {"Intra": intra, "freq": freq, 'Test': test, "Type": types, "Category": category}
df = pd.DataFrame(data)
df.to_csv(os.path.join(path, "existence.csv"))
plot_intra_perc_violin(df, os.path.join(path, 'max_freq_violin.svg'), "Test", "Intra", "Type")

dispa= np.repeat(x,z)
nodispa= np.repeat(x2,z2)
heada= np.repeat(x3,z3)
paleoa= np.repeat(x4,z4)

abs=[dispa, nodispa, heada, paleoa]
freq=[]
for i in abs:
    freq.extend(i)

test=[]
types=[]
for i in dispa:
    test.append('Disp')
    types.append('Abs')
for i in nodispa:
    test.append('No Disp')
    types.append('Abs')
for i in heada:
    test.append('Head')
    types.append('Abs')
for i in paleoa:
    test.append('Paleo')
    types.append('Abs')

df=pd.DataFrame()
df['Intra']=freq
df['Test']=test
df.to_csv(os.path.join(path,"freq_abs.csv"))

subprocess.call ("/usr/bin/Rscript --vanilla mat_plots.R", shell=True)
