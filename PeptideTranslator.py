#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 19:12:09 2019

@author: Colin Tang

Script to translate unique peptide library to nucleotide sequence and 
tests sequences for levenshtein distance.

inputs are .txt files. AACodonFile is tab delimited with codons comma delimited
"""
    
import csv
import random
from Levenshtein import distance
import datetime
    
def PeptideTranslator(peptideLibFile, AACodonFile, levenshteinLimit = 3):
    #builds peptide library list
    with open(peptideLibFile) as f:
        peptideLib = f.read().splitlines()
        f.close()
    
    #builds AA codon dictionary
    with open(AACodonFile) as a: 
        AAcodonList = a.read().splitlines()
        a.close()
    
    AAcodonDict = {}
    for i in AAcodonList:
        AArow = i.split("\t")
        AAcodonDict[AArow[0]] = AArow[1].split(",")
        
    #performs first nucleotide translation
    peptideDict = {}
    aaList = list(peptideLib[0])
    translatedSeq = list()
    #loop to translate peptide sequence to nucleotide codon randomly
    for i in aaList:
        translatedSeq.append(random.choice(AAcodonDict[i]))
    
    peptideDict[peptideLib[0]] = "".join(translatedSeq)
    peptideNum = 1   #counter for next peptide translation
    del peptideLib[0]   #removes the first peptide from list
    
    #loop to translate the rest of the peptides and perform the levenshtein tests
    for i in peptideLib:
        peptideNum += 1   #increases the counter for current peptide
        aaList = list(i)
        #loop to translate each amino acid
        while len(peptideDict) < peptideNum:
            translatedSeq = list()
            for j in aaList:
                translatedSeq.append(random.choice(AAcodonDict[j]))
            
            translatedSeq = "".join(translatedSeq)
            peptideDictSeq = list(peptideDict.values())
            #loop to get levenshtein score for every nucleotide sequence in the dictionary
            levenshteinScores = list()
            for j in peptideDictSeq:
                levenshteinScores.append(distance(translatedSeq, j))
                
            levenshteinScores.sort()
            if levenshteinScores[0] > levenshteinLimit:
                peptideDict[i] = translatedSeq
                if peptideNum % 100 == 0:
                    print("Translated " + str(peptideNum) + " peptide sequences at " + str(datetime.datetime.now()))
        
    return peptideDict
        
        
    
#function to write dictionary x to file

def dictWriter(outputDict, outputName = "dictname.csv"):
    with open(outputName, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in outputDict.items():
           writer.writerow([key, value])

    
    
    
    
