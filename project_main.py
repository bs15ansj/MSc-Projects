#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 14:05:57 2020

@author: bs15ansj
"""

import consensus
import argparse


def Main():
    parser = argparse.ArgumentParser()
    
    # Parser for the main commands arguments
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-b', '--blast', help="Run blast",
                        action='store_true')
    group.add_argument('-a', '--align', action='store_true', 
                        help = "Align sequences from xml or fasta file")
    group.add_argument('-c', '--consensus', action='store_true', 
                        help = "Caclulate the consensus sequence from an alignment")
    group.add_argument('-ba', '--blast_align', action='store_true', 
                        help = "Run blast then align using clustal")
    group.add_argument('-bc', '--blast_consensus', action='store_true', 
                        help = "Run blast then align and caclulate the consensus sequence")
    group.add_argument('-pc', '--plot_consensus', action='store_true', 
                        help = "Plot the consensus sequence from an alignment")

    # File name arguments
    # Output file names can be specified using the following arguments. If none
    # are provided then the default will be the name of the input file with 
    # an added extension e.g. input -> input.blast_results.xml
    parser.add_argument('-sf', '--sequence_file', type=str, 
                        help="Name of sequence file.", metavar='')
    parser.add_argument('-bf', '--blast_file', type=str,
                        help="Name of blast file. Must have .xml extension.",
                        metavar='')
    parser.add_argument('-af', '--alignment_file', type=str,
                        help="Name of alignment file. Must have .aln extension.",
                        metavar='')
    parser.add_argument('-ff', '--fasta_file', type=str,
                        help="Name of blast output fasta file." \
                             " Must have .fasta extension.",
                        metavar='')
    parser.add_argument('-cf', '--consensus_file', type=str,
                        help="File to write consensus sequence to.",
                        metavar='')
    parser.add_argument('-pf', '--plot_file', type=str,
                        help="File to write plot to.",
                        metavar='')    
    # Options
    parser.add_argument('-mh', '--max_hits', type=int, help='Maximum number '\
                        'of hit sequences for blast to return. Default = 100', 
                        default=100)
        
    args = parser.parse_args()
    
    
    # Run blast if --blast, --blast_align, or --blast_consensus are used. 
    if (args.blast or args.blast_align or args.blast_consensus) and (args.sequence_file is None):
        print("\nError!:\n\tBlast requires a sequence: use --sequence_file")
        return
    if (args.blast or args.blast_align or args.blast_consensus) and (args.blast_file is None):
        args.blast_file = args.sequence_file.split('.')[0]+'.blast_results.xml'
    if args.blast or args.blast_align or args.blast_consensus:
        consensus.Run_blast(infile=args.sequence_file, 
                            outfile=args.blast_file, 
                            max_hits=args.max_hits)

    # Align
    if args.align:
        if args.blast_file is None:
            if args.fasta_file is None:
                print("\nError!:\n\tAlign requires sequences: use --blast_file or "\
                  "--fasta_file")
                return
            else:
                if args.alignment_file is None:
                    args.alignment_file = args.fasta_file.split('.')[0]+'.clustal_alignment.aln'
                consensus.Clustal_alignment(fastafile=args.fasta_file, 
                                            alnfile=args.alignment_file)
                return
        else:
            if args.alignment_file is None:
                    args.alignment_file = args.blast_file.split('.')[0]+'.clustal_alignment.aln'
            consensus.Clustal_alignment(xmlfile=args.blast_file, 
                                        alnfile=args.alignment_file)
            return
    
    # Blast align and blast consensus
    if (args.blast_align or args.blast_consensus) and (args.alignment_file is None):
        args.alignment_file = args.sequence_file.split('.')[0]+'.clustal_alignment.aln'
    if (args.blast_align or args.blast_consensus) and (args.fasta_file is None):
        consensus.Clustal_alignment(xmlfile=args.blast_file, 
                                    alnfile=args.alignment_file)
    elif args.blast_align or args.blast_consensus:
        consensus.Clustal_alignment(xmlfile=args.blast_file, 
                                    fastafile=args.fasta_file,
                                    alnfile=args.alignment_file)
    
    # Calculate consensus if --blast_consensus or --consensus are used.
    if args.blast_consensus and (args.consensus_file is None):
        args.consensus_file = args.sequence_file.split('.')[0]+'.consensus.txt'
    if args.consensus and (args.alignment_file is None):
        print("\nError!:\n\tConsensus calculation requires multiple sequence"\
              " alignment: use --alignment_file")
        return
    if args.consensus and (args.consensus_file is None):
        args.consensus_file = args.alignment_file.split('.')[0]+'.consensus.txt'
    if args.blast_consensus or args.consensus:
        consensus.Consensus_seq(alnfile=args.alignment_file,
                                confile=args.consensus_file)
        
    # Plot consensus if --plot_consensus is used
    if args.plot_consensus:
        if args.alignment_file is None:
            print("\nError!:\n\tConsensus plotting requires multiple sequence"\
                  " alignment: use --alignment_file")
            return
        elif args.plot_file is None:
            args.plot_file = args.alignment_file.split('.')[0]+'.png'
            consensus.Plot_consensus(alnfile=args.alignment_file,
                                     plotfile=args.plot_file)
        else:
            consensus.Plot_consensus(alnfile=args.alignment_file,
                                     plotfile=args.plot_file)        
    return
        
    

if __name__ == '__main__':
    Main()