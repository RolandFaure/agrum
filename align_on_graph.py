#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 17:45:51 2022

@author: rfaure
"""

import segment

#recursive function to find a path between two contigs
#input: two contigs, the graph, and a range of acceptable distance between the two contigs
#output: the path linking the two contig if it is the only one that has an acceptable length. If there are no paths or several, return -1
def find_path(segments, contig1, end1, contig2, end2, min_distance, max_distance, length_now = -1, path_now = "", already_seen_contigs={}):
    
    #print("Entering: ", contig1.names, " ", end1, " ", contig2.names, " ", end2, " ", min_distance, " ", max_distance, " ", length_now, " ", path_now)
    if (contig1, end1) in already_seen_contigs:
        if already_seen_contigs[(contig1, end1)] == True: #another victorious path exists from here
            # print("Contig: ", contig1.names, " ", end1, " ", contig2.names, " ", end2, " ", min_distance, " ", max_distance, " ", length_now, " ", path_now)

            return ["no", "paths"] #of length 2, making it impossible to find a single path
        else : #there is no path to victory through there and we know it
            return []
        
    allpaths = []
    if length_now > max_distance:
        return []
    if contig1 == contig2 and end1 == 1-end2 and length_now >= min_distance : #we found a path !
        allpaths += [path_now]
       
    if length_now == -1:
        length_now = 0
    else:
        length_now+=contig1.length
        
    for n, neighbor in enumerate(contig1.links[end1]):
        orientation = "+,"
        if contig1.otherEndOfLinks[end1][n] == 1 :
            orientation = "-,"
        allpaths += find_path(segments, neighbor, 1-contig1.otherEndOfLinks[end1][n], contig2, end2, min_distance, max_distance, length_now, path_now+neighbor.names[0]+orientation, already_seen_contigs)
    
    already_seen_contigs[(contig1, end1)] = (len(allpaths)>0)

    return allpaths
    
    

#input: segments, names (gfa graph), alignments (fasta file of reads aligned on graph), out (out file)
def align_on_graph(segments, names, file_alignments, outfile):
        
    alignments = open(file_alignments, "r")
    out = open(outfile, "w")
    current_read = ""
    paf_fields = []
    
    for line in alignments:
        
        ls = line.split("\t")
        if ls[0] == current_read:
            
            everythingallright = True
            contig = 0
        
            if len(ls) > 11 and ls[5] in names :
                contig = ls[5]              
            else :
                everythingallright = False
                
            if everythingallright and len(ls) > 11:
                paf_fields += [(int(ls[2]), int(ls[3]), ls[4], contig, int(ls[7]), int(ls[8]), int(ls[9]), int(ls[10]), int(ls[11]))]  # (start on the read, end on the read, orientation, contig, start on contig, end on contig, #of residue matches, alignment block length, mapping quality)
        
        else: #we change read, write the alignment in the augmented paf
        
            if len(paf_fields)>= 1:
                
                paf_fields.sort(key= lambda x: x[0])
                
                path = {}
                path["start_on_read"] = str(paf_fields[0][0])
                path["start_on_first_contig"] = str(paf_fields[0][4])
                path["path"] = str(paf_fields[0][3])+paf_fields[0][2]+","
                for contigaln in range(0, len(paf_fields)-1):
                    
                    # path["path"] += str(paf_fields[contigaln][3])+paf_fields[contigaln][2]+","
                    
                    min_distance = float(paf_fields[contigaln+1][0] - paf_fields[contigaln][1])/2 - 300
                    max_distance = max( float(paf_fields[contigaln+1][0] - paf_fields[contigaln][1])*2, 1000 )
                    
                    end1 = 1
                    if paf_fields[contigaln][2] == "-":
                        end1 = 0
                    end2 = 0
                    if paf_fields[contigaln+1][2] == "-":
                        end2 = 1
                    
                    # if paf_fields[contigaln][3] == "193361":
                    #print("between contig ", paf_fields[contigaln][3], " and ", paf_fields[contigaln+1][3], " input: ", end1, " ", end2, " ", min_distance, " ", max_distance)
                    contigs_between = find_path(segments, segments[names[paf_fields[contigaln][3]]], end1, segments[names[paf_fields[contigaln+1][3]]], end2, min_distance, max_distance, already_seen_contigs={})
                    #print("path between", paf_fields[contigaln][3], " and ", paf_fields[contigaln+1][3], ", found these paths: ", contigs_between)
            
                    if len(contigs_between) == 1:
                        path["path"] += contigs_between[0]
                    else : #means the two contigs cannot be linked
                        alignmentLine = current_read+"\t"+path["start_on_read"]+"\t"+str(paf_fields[contigaln][1])+"\t"+path["start_on_first_contig"]+"\t"+str(paf_fields[contigaln][5])+"\t"+path["path"].rstrip(",")+"\n"
                        out.write(alignmentLine)
                        path["start_on_read"] = str(paf_fields[contigaln+1][0])
                        path["start_on_first_contig"] = str(paf_fields[contigaln+1][4])
                        path["path"] = str(paf_fields[contigaln+1][3])+paf_fields[contigaln+1][2]+","                       
                
                alignmentLine = current_read+"\t"+path["start_on_read"]+"\t"+str(paf_fields[-1][1])+"\t"+path["start_on_first_contig"]+"\t"+str(paf_fields[-1][5])+"\t"+path["path"].rstrip(",")+"\n"
                out.write(alignmentLine)
        
            #now start with a new read
            paf_fields = []
            current_read = ls[0]
            
            ls = line.split("\t")
            if ls[0] == current_read:
                
                everythingallright = True
                contig = 0
        
            if len(ls) > 11 and ls[5] in names :
                contig = ls[5]                
            else :
                everythingallright = False
                
            if everythingallright and len(ls) > 11:
                paf_fields += [(int(ls[2]), int(ls[3]), ls[4], contig, int(ls[7]), int(ls[8]), int(ls[9]), int(ls[10]), int(ls[11]))]  # (start on the read, end on the read, orientation, contig, start on contig, end on contig, #of residue matches, alignment block length, mapping quality)
        
    #last read
    for contigaln in paf_fields:
        alignmentLine += str(contigaln[3]) + str(contigaln[2])+","
    alignmentLine.rstrip(",")
    alignmentLine += "\t"