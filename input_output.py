#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:59:45 2020

File basically dedicated to small functions involving reading and writing files
"""
import numpy as np
from scipy import sparse #to handle interactionMatrix, which should be sparse
import time #to inform the user on what the programm is doing on a regular basis
import os.path #to check the existence of files
import pickle #for writing files and reading them
import re #to find all numbers in a mixed number/letters string (such as 31M1D4M), to split on several characters (<> in longReads_interactionMatrix)
import shutil #to remove directories
import sys #to exit when there is an error and to set recursion limit


from segment import Segment
from segment import compute_copiesNumber
from segment import delete_links_present_twice



# Return a list in which each element contains a list of linked contigs (accroding to GFA). There is one list for each end of the contig
# Also returns the list of the contig's names
def load_gfa(file):

    #print('Loading contigs')
    gfa_read = open(file, "r")

    segments = []
    
    index = 0
    names = {} # names is a dictionary that associates the name of each contig in the gfa with an index (which will correspond later to the one in interactionMatrix and copiesnumber)
    
    for line in gfa_read:
        if line[0] == "S":
            l = line.strip('\n').split("\t")
            cov = 0
            
            for element in l :
                if 'dp' in element[:2] or 'DP' in element[:2] or 'rd' in element[:2] :
                    try :
                       cov = float(element.split(":")[-1])
                    except:
                        pass
                        
                elif 'RC' in element[:2] or 'KC' in element[:2] :
                    try :
                       cov = float(element.split(":")[-1])/len(l[2])
                    except:
                        pass
            
            s = Segment([l[1]], [1], [len(l[2])], readCoverage = [cov])
            segments.append(s)
            names[s.names[0]] = index #now this contig (identified by its name) is attached to index
            index += 1
            

    #print('Loading links')
    gfa_read = open(file, "r")
        
    cov = 1
    for line in gfa_read:
        if line[0] == "L":

            l = line.strip('\n').split("\t")
            
            segments[names[l[1]]].add_link_from_GFA(line, names, segments, 0)
            segments[names[l[3]]].add_link_from_GFA(line, names, segments, 1)

    gfa_read.close()
    
    delete_links_present_twice(segments)

    return segments, names

