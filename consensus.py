
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 20:51:30 2020
@author: bs15ansj
This is a small script to help you get used to using biopython. In order to
perform the clustalw multiple sequence alignment, make sure you have clustal
installed. You can do this using the following commands:
    
    sudo apt-get update
    sudo apt-get install clustalw
    
Once installed, make sure that you close the terminal and reopen before 
running the script. 
"""

from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

aas =['V','I','L','E','Q', \
'D','N','H','W','F','Y', \
'R','K','S','T','M','A', \
'G','P','C','-']

def Run_blast(local=False, infile=None, outfile=None, 
              max_hits=100):

    # Read the input file and get the sequence
    print('\nBlasting sequence from '+infile+'...')
    print('Maximum number of hits set to '+str(max_hits))
    with open(infile, "r") as f:
        in_file = f.read().strip()

    if local == False:
        
        # Blast the sequence over the internet (you may have to wait) and 
        # save the results in a qblast object
        result_handle = NCBIWWW.qblast("blastp", "nr", 
                                       in_file, hitlist_size=max_hits)
        
        
        # Open the output file then read the results object and write to 
        # the output file
        print('\tDone: writing to '+outfile)
        
        with open(outfile, "w") as out_handle:
            out_handle.write(result_handle.read())
            
    
    if local == True:
        
        cmd = NcbiblastpCommandline(query=in_file, db='nr', outfmt=5, 
                                    out=outfile)
        cmd
        
def xml2fasta(infile=None, outfile=None):
    
    print('\nConverting '+infile+' to fasta format, removing duplicates...')

    
    # Load the blast output file
    blast_qresult = SearchIO.read(infile, 
                                  "blast-xml")

    # Iterate through ids and sequences and add them to lists. Sequences are
    # only added to the list if they are not already in the list. This is done
    # because sometimes there will be duplicates in the blast output, which 
    # will produce an error when aligning with clustal 
    
    ids = []
    sequences = []
    for hsp in blast_qresult.hsps:
        if hsp.hit.id not in ids:
            ids.append(str(hsp.hit.id))
            sequences.append(str(hsp.hit.seq))
            
    # Open the sequences output file then for each high-scoring pair
    # in the blast results, write the hit ID (proceeded by a ">" for 
    # fasta format), followed by the hit sequence on the next line    
    with open(outfile, "w") as f:
        for i, s in zip(ids, sequences):
            f.write('> '+i+'\n')
            f.write(s+'\n')
    
    print('\tDone: writing to '+outfile)

def Clustal_alignment(xmlfile=None, fastafile=None, alnfile=None):
    
    if fastafile is None:
        fastafile = xmlfile.replace('.xml','.fasta')
        xml2fasta(infile=xmlfile, outfile=fastafile)
    
    # Run the command line version of clustal using the sequences.fasta
    # file and output to a clustal format alignment file
    print('\nAligning '+fastafile+' with clustal...')
    clustalo_exe = r"D:\Affimer_project\clustal-omega-1.2.2-win64\clustalo.exe"
    cmd = ClustalwCommandline(clustalo_exe, 
                              infile=fastafile, 
                              outfile=alnfile)
    cmd()
    print('\tDone: writing to '+alnfile)

    
def MSA_matrix(alnf):
    
    # Read the clustal alignment
    align = AlignIO.read(alnf,'clustal')
    
    # Read sequences, convert each to list and add to matrix
    matrix = []
    for record in align:
        matrix.append(list(record.seq))
    matrix = np.array(matrix)
    return matrix     
    
def frequency(alnf):

    # Generate the MSA matrix and create a new empty frequency matrix
    
    M = MSA_matrix(alnf)
    F = []
    
    # Set the number of sequences equal to the number of rows in the matrix
    # (using .shape on a matrix returns a two integer tuple, the first number 
    # is the number of rows, the second is the number of colums)
    
    no_sequences = M.shape[0]
    
    # Enumerate through the transpose of the MSA matrix. 
    # The transpose of a matrix flips the matrix along the diaganol axis e.g.
    # 
    #  A B C               A D G    
    #  D E F  transpose => B E H    
    #  G H I               C F I
    # 
    # This is because when we enumerate, we go through the matrix row by row. 
    # We want to enumerate through the columns so we transpose the matrix, 
    # effectively switching the rows and colums. 
    
    for index, i in enumerate(M.T):
        
        # Create an empty list that will store the frequency of each amino
        # acid for position i
        
        f = []
        
        # Iterate over the amino acids and count their frequency in position i
        # and add this to the list
        
        for aa in aas:
            x = np.sum(np.char.count(i, aa)) / no_sequences
            f.append(x)
        
        # After creating the list of frequencies, f, add this to the frequency
        # matrix, F, for each position i
            
        F.append(f)
    
    # Return the frequency matrix as an array  
        
    return np.array(F)
            
def Consensus_seq(alnf,alnfile=None, confile=None):
    
    print('\nCalculating consensus sequence from '+alnfile+'...')
    
    
    # Generate a frequency matrix
    F = frequency(alnf)
    
    # Create empty lists for the consensus sequence, and the frequency of 
    # consensus residues
    con_seq = []
    con_freq = []
    
    # Iterate over rows of the matrix (remember, this has been transposed so 
    # this is equivalent to iterating over columns of the MSA)
    for i in F:
        
        # Find column where the frequency is maximum. This column index
        # corresponds to the list index for a particular amino acid in aas at 
        # the top of this script. 
        m = np.where(i == np.amax(i))
        
        # Use the column index in m to find the amino acid letter in aas, then
        # add this letter to the consensus sequence
        con_seq.append(aas[int(m[0][0])])
        
        # Add the max frequency to the consensus frequency list
        con_freq.append(np.amax(i))
    
    if confile:
        with open(confile, 'w+') as f:
            for aa, i in zip(con_seq, con_freq):
                f.write(aa+' '+str(i)+'\n')
        f.close()
            
        print('\tDone: writing to '+confile)

    return con_seq, con_freq

def Plot_consensus(alnf,x,alnfile=None, plotfile=None):
    
    # Generate the consensus frequency list (consensus_seq() returns two
    # variable, so we need to create a variable name for both even though
    # we will only need con_freq)
    con_seq, con_freq = Consensus_seq(alnf,alnfile=alnfile)
    print('\tDone')
    print('\nPlotting consensus sequence from '+alnfile+'...')
    plt.plot(range(1, len(con_freq)+1), con_freq, label='Hits: '+str(x))
    plt.xlim(1, len(con_freq))
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.legend(bbox_to_anchor=(1.10,1))
    plt.savefig(plotfile, dpi = 300)
    print('\tDone: writing to '+plotfile)
    f = open('Consensus Freq.txt','a+')
    f.write('Number of hits : ' + str(x) +'\n')
    f.write('Frequency: '+str(con_freq)+'\n')
    f.write("Sequence: "+str(con_seq)+'\n')

def Blast_a_Bunch(x,itter=5):
    for i in range(itter):
        alnf = str(x)+'_align.aln'
        Run_blast(infile='input.fasta',outfile='BLAST'+str(x)+'.xml',max_hits=x)
        xml2fasta(infile='BLAST'+str(x)+'.xml',outfile='Blast_'+str(x)+'.fasta')
        Clustal_alignment(fastafile='Blast_'+str(x)+'.fasta',alnfile=alnf)
        Plot_consensus(alnf,x,alnfile=alnf,plotfile='total_graphs.png')
        x = x+100
