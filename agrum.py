#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:02:00 2022

@author: rfaure
"""

import argparse
from input_output import load_gfa
from align_on_graph import align_on_graph
import sys
import os
SORT = True

def printRed(skk): print("\033[91m {}\033[00m" .format(skk))
def printGreen(skk): print("\033[92m {}\033[00m" .format(skk))

def check_dependancies():
    
    code = os.system("minimap2 -h > trash.txt 2> trash.txt")
    if code != 0 :
        print("WARNING: minimap2 command does not answer. Install minimap2 (on the PATH) or use option -a")
    
    code = os.system("awk -h > trash.txt 2> trash.txt")
    if code != 0 :
        print("WARNING: awk command does not answer. Install awk or use option -a")
        
    code = os.system("sort --help > trash.txt 2> trash.txt")
    if code != 0 :
        print("WARNING: sort command does not answer. Install sort or use option -a with a file SORTED BY READ NAMES")
        SORT = False

def parse_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g",
        "--graph",
        required=True,
        help="""Assembly graph (GFA format)""",
    )
    parser.add_argument(
        "-x",
        "--preset",
        required=True,
        help="""Presets for minimap2. Accepted values are map-ont (Nanopore), map-pb (Pacbio CLR) and map-hifi (HiFi)""",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="""Output file (.tsv format) """,
    )
    parser.add_argument(
        "-r",
        "--reads",
        required=False,
        default = "Empty",
        help="""Long reads""",
    )
    parser.add_argument(
        "-a",
        "--alignments",
        required=False,
        default = "Empty",
        help="""Reads aligned on the fasta contigs (PAF format)""",
    )
    
    return parser.parse_args()

def main():
    
    print(   """
            .....
        _d^^^^^^^^^b_
     .d''           ``b.
   .p'    Welcome in   `q.
  .d'  Agrum (Aligning  `b.
 .d'   on Graphs Using   `b.
 ::       Minimap2)       ::
 ::  ...................  ::
 ::                       ::
 `p.     A sofware by    .q'
  `p.    Roland Faure   .q'
   `b.    12.10.2022   .d'
     `q..          ..p'
        ^q........p^
            ''''
            """)
    
    check_dependancies()
    
    args = parse_args()
    gfaFile = args.graph
    alignments = args.alignments
    readsFile = args.reads
    preset = args.preset
    outfile = args.out
        
    if not (preset == "map-ont" or preset == "map-pb" or preset == "map-hifi" ):
        print("ERROR: please choose a valid preset (-x option), among map-ont, map-pb or map-hifi")
    
    if readsFile == "Empty" and alignments == "Empty":
        
        print("ERROR: provide either long reads or alignments")
        sys.exit()
        
    elif alignments == "Empty":
        
        printRed("====== STAGE 1 ======\n")
        printGreen("--- Step 1 - convert GFA file to fasta file ---\n ")
        
        fastaname = gfaFile.rstrip(".gfa")+ ".fa "
        command = "awk '/^S/{print \">\"$2\"\\n\"$3}' " + gfaFile + " > " + fastaname
        print("To convert, Agrum will use awk command:\n\t ", command)
        os.system(command)
        
        print("GFA file converted in the fasta file ", fastaname)
        
        printGreen("\n--- Step 2 - using minimap2 to map the reads on the fasta ---\n")
        
        command = "minimap2 --secondary=no -x " + preset + " " + fastaname + " " + readsFile + " > alignments.paf 2> trash.txt"
        print("To map, Agrum will use command:\n\t", command)
        os.system(command)
        alignments = "alignments.paf"
         
        print("Alignments are stored in file alignments.paf. If you reuse agrum on the same dataset, directly feed it this file through the -a option")
        
    else:
         printRed("====== STAGE 1 ======\n\n")
         print("Skipped because the alignment file is already provided through the -a option")
        
    if SORT: #check if alignment file is sorted
    
        printRed("\n\n====== STAGE 2 ======\n\n")
        print("Agrum needs the alignment file to be sorted by read name. Making sure it is.")
        
        code = os.system("sort -c -n -k 1,1 -s " + alignments) #check if sorted
        
        if code != 0: #if not sorted, then sort
            command = "sort " + alignments + " -n -k 1,1 -s -S 5G > tmp.paf"
            print("The alignment file is not sorted ! Sorting it with command:\n\t", command, " & mv tmp.paf ", alignments)
            os.system(command)
            os.system("mv tmp.paf " + alignments)
            
        else:
            print("Alignment file is already sorted! Proceeding to stage 3.")
    else:
        printRed("\n\n====== STAGE 2 ======\n\n")
        print("Agrum needs the alignment file to be sorted by read name. However, sort command does not seem to be installed on this computer. This should be ok, but you can check manually if the .paf is indeed sorted by read name")
    
    printRed("\n\n====== STAGE 3 ======\n\n")
    
    print("In this stage, I will try to infer alignent paths from alignments on the contigs. There are essentially two possibilities: ")
    
    print("CASE 1: ")
    print("""
     contig 1                                            contig 4
=================>.                               .===================>
                   '.	      contig 3	        .'
                     <==========================
    contig 2       .'                           '.       contig 5
=================>'                               '<===================

          """)
    print("Read x mapped on contig 1 and contig 4 -> Agrum infers that x maps on the path {contig 1+, contig 3-, contig 4+}\n")
    print("CASE 2: ")
    print("""
                               contig 2
                     .<=========================.
     contig 1      .'                            '.       contig 4
=================>:                                :===================>
                   '.          contig 3          .'
                     '=========================>'

          """)
    print("Read x mapped on contig 1 and contig 4 -> Agrum cannot infer the path (it could map either on contig 2 or contig 3) :/\nFor that case, you need a real mapper on graph, such as minichain, GraphAligner, Spaligner...")
    
    print("\nOutput will be written in ", outfile, " as a tsv file. It is a tab-separated file in which the fields are: ")
    
    print("""
     _______________________________________________________________
    | Field |                         Value                         |
    |-------|-------------------------------------------------------|
    |   1   |                    Name of the read                   |
    |   2   |      First position of the alignment on the read      |
    |   3   |      Last position of the alignment on the read       |
    |   4   |  First position of the alignment on the first contig  |
    |   5   |   Last position of the alignment on the last contig   |
    |   6   |   Path of the alignement (e.g. ctg_1+,ctg_3-,ctg_4+)  |
     ------- -------------------------------------------------------
          """)
    print("Currently infering paths...")
    
    segments, names = load_gfa(gfaFile)
    align_on_graph(segments, names, alignments, outfile)
    
    os.system("rm trash.txt")
    
    print("Done! Hope you had a great time running Agrum! If you've had a problem or would like additional features on Agrum, drop me an issue on the GitHub")
    

if __name__ == "__main__":
    main()