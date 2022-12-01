# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 08:55:37 2022

@author: gm
"""

import pandas as pd
import pyreadstat
import numpy as np
import matplotlib.pyplot as plt


class Gene():   
    def __init__(self):
        
        self.neo=[]
        self.type=[]
        self.name=""
        self.count=0
        self.index=0
        self.loaded=False
                      
    # name of the gene    
    def setName(self, name):
        self.name=name
        
    # Mutation based on neoclotide such as C->G
    def setNeo(self, neo):
        self.neo.append(neo)
    
    # mutation type including Missense or Silent    
    def setType(self, type):
        self.type.append(type)

    # counts of each gene
    def setCount(self, count):   
        self.count=count

    # index of gene in the list ,
    def setIndex(self, index):   
        self.index=index
        self.loaded=True


# =============================================================================
#          # Finding the mutations; For example input: 'GA>GG' , output: 'AG'  
# =============================================================================

def find_neo(s):
    
    s1=s[0:2]
    s2=s[3:5]
    
    overlap=set(s1).intersection(s2)
    
    if (len (overlap)!=0):
        
        a=list(overlap)[0]
        s1_index=s1.index(a)
        s2_index=s2.index(a)
        
        if (s1_index==0 and s2_index==0):
            mut_index=1
        else:
            mut_index=0
            
        mut=s1[mut_index]+s2[mut_index]
        
    else:
        
        mut=s1[0]+s2[0]
        
        
    return mut


# =============================================================================
#          # Finding the number of mutations of each neo   
# ============================================================================= 
def count_mutationNeo(first_g_neo_mut):
            # number of each mutation type in first gene
            
    neos_count=[]
    all_neos=['CT','TA','AG','CA','GT','GA','TG','CG','AC','TC','GC','AT'] 
    
    for i in range(len(first_g_neo_mut)):
        for k in range (len(all_neos)):
            neos_count.append(first_g_neo_mut.count(all_neos[k]))
        
   
    result={'CT':neos_count[0], 'TA':neos_count[1],'AG':neos_count[2],'CA':neos_count[3],'GT':neos_count[4],'GA':neos_count[5],
            'TG':neos_count[6],'CG':neos_count[7],'AC':neos_count[8],'TC':neos_count[9],'GC':neos_count[10],'AT':neos_count[11],}
    return result


# =============================================================================
#          # Finding the number of mutations of each type   
# ============================================================================= 
def count_mutationType(gene):
            # number of each mutation type in first gene
    missense_count=gene.type.count('Missense_Mutation')
    nonsens_count=gene.type.count('Nonsense_Mutation')
    rna_count=gene.type.count('RNA')
    silent_count=gene.type.count('Silent')
    splice_count=gene.type.count('Splice_Site')
    nonstop_count=gene.type.count('Nonstop_Mutation')
    
    result={'Missense_Mutation':missense_count, 'Nonsense_Mutation':nonsens_count,'RNA':rna_count, 'Silent':silent_count,'Splice_Site': splice_count, 'Nonstop_Mutation':nonstop_count }
    return result

# =============================================================================
#          # calculate similarity w.r.t type   
# =============================================================================    
def  similarity_type(first_g, sec_g):  

    
     all_types=['Missense_Mutation','Nonsense_Mutation','RNA','Silent','Splice_Site','Nonstop_Mutation' ] 

     type_intersec=[]
     type_intersec=list(set(first_g.type).intersection(sec_g.type))
      
     # w subset: the intersection classes
     w_sub= type_intersec
     cw_sub= [i for i in all_types if i not in w_sub]
     
     #loop through the w subset             
     w_len=len(w_sub) 
     cw_len=len(cw_sub) 

# =============================================================================
#  simlarity between first and second: sigma (first ,second)
# =============================================================================
     # calculating p(w)  

     if w_len!=0:     # if there are some mutation types in common between first and second genes
                 
      # finding how many mutations of the same type of from the common mutations we have in the first gene  
      for t in range (len(first_g.type)):                                    
            count_mut_dic_first=count_mutationType(first_g)

# =============================================================================
# # Finding the sum pf probabilities for first gene (p(w) in (first,sec) order)                
# =============================================================================
            
            count_mut_dic_first=count_mutationType(first_g)
            
            prob_w=0                # probability for all mutations for first gene in (first,sec) order
            for j in range(w_len):
                current_type=w_sub[j]
                if (current_type=="Missense_Mutation"):
                    prob_w=prob_w + (count_mut_dic_first['Missense_Mutation']/first_g.count)
                    
                if (current_type=="Nonsense_Mutation"):
                        prob_w=prob_w + (count_mut_dic_first['Nonsense_Mutation']/first_g.count)
                        
                if (current_type=="RNA"):
                              prob_w=prob_w + (count_mut_dic_first['RNA']/first_g.count)
                              
                if (current_type=="Silent"):
                              prob_w=prob_w + (count_mut_dic_first['Silent']/first_g.count)
                                    
                if (current_type=="Splice_Site"):
                              prob_w=prob_w + (count_mut_dic_first['Splice_Site']/first_g.count) 
                              
                if (current_type=="Nonstop_Mutation"):
                             prob_w=prob_w + (count_mut_dic_first['Nonstop_Mutation']/first_g.count) 
                              
                              
# =============================================================================
#  Finding the sum pf probabilities for second gene (p(cw) in (first,sec) order)  
# =============================================================================
           
            count_mut_dic_sec=count_mutationType(sec_g)        
            prob_cw=0                # probability for all mutations for first gene in (first,sec) order
            for j in range(cw_len):
                current_type=cw_sub[j]
                if (current_type=="Missense_Mutation"):
                    prob_cw=prob_cw + (count_mut_dic_sec['Missense_Mutation']/sec_g.count)
                    
                if (current_type=="Nonsense_Mutation"):
                        prob_cw=prob_cw + (count_mut_dic_sec['Nonsense_Mutation']/sec_g.count)
                        
                if (current_type=="RNA"):
                              prob_cw=prob_cw + (count_mut_dic_sec['RNA']/sec_g.count)
                              
                if (current_type=="Silent"):
                              prob_cw=prob_cw + (count_mut_dic_sec['Silent']/sec_g.count)
                                    
                if (current_type=="Splice_Site"):
                              prob_cw=prob_cw + (count_mut_dic_sec['Splice_Site']/sec_g.count)    
                                
                if (current_type=="Nonstop_Mutation"):
                             prob_cw=prob_cw + (count_mut_dic_sec["Nonstop_Mutation"]/sec_g.count)               

# =============================================================================
# # Overall  similarity based on type between first and second gene 
# =============================================================================
            sigma_type= prob_w +prob_cw - 1         # sigma= p(w) +p(cw) - 1 
            
     else:
         sigma_type=1
       
     return sigma_type  


# =============================================================================
#          # calculate similarity w.r.t neo   
# =============================================================================    

def  similarity_neo(first_g, sec_g):  
     
     all_neos=['CT','TA','AG','CA','GT','GA','TG','CG','AC','TC','GC','AT']  
     type_intersec=[]
     first_g_neo_mut=[]
     sec_g_neo_mut=[]
    
     for i in range (len(first_g.neo)):
         first_g_neo_mut.append(find_neo(first_g.neo[i]))
         
     for i in range (len(sec_g.neo)):
         sec_g_neo_mut.append(find_neo(sec_g.neo[i]))
    
     type_intersec=list(set(first_g_neo_mut).intersection(sec_g_neo_mut))
        

     # w subset: the intersection classes
     w_sub= type_intersec
     cw_sub= [i for i in all_neos if i not in w_sub]
     
     #loop through the w subset             
     w_len=len(w_sub) 
     cw_len=len(cw_sub) 

# =============================================================================
#  simlarity between first and second: sigma (first ,second)
# =============================================================================
     # calculating p(w)  

     if w_len!=0:     # if there are some mutation types in common between first and second genes
                 
      # finding how many mutations of the same type of from the common mutations we have in the first gene  
         #for t in range (len(first_g_neo_mut)):                                    
        count_mut_dic_first=count_mutationNeo(first_g_neo_mut)

# =============================================================================
# # Finding the sum pf probabilities for first gene (p(w) in (first,sec) order)                
# =============================================================================

        prob_w=0                # probability for all mutations for first gene in (first,sec) order
        for j in range(w_len):
            current_type=w_sub[j]
            if (current_type=="CT"):
                prob_w=prob_w + (count_mut_dic_first['CT']/first_g.count)
                
            if (current_type=="TA"):
                    prob_w=prob_w + (count_mut_dic_first['TA']/first_g.count)
                    
            if (current_type=="AG"):
                          prob_w=prob_w + (count_mut_dic_first['AG']/first_g.count)
                          
            if (current_type=="CA"):
                          prob_w=prob_w + (count_mut_dic_first['CA']/first_g.count)
                                
            if (current_type=="GT"):
                          prob_w=prob_w + (count_mut_dic_first['GT']/first_g.count) 
                          
            if (current_type=="GA"):
                          prob_w=prob_w + (count_mut_dic_first['GA']/first_g.count)
                    
            if (current_type=="TG"):
                          prob_w=prob_w + (count_mut_dic_first['TG']/first_g.count)
                          
            if (current_type=="CG"):
                          prob_w=prob_w + (count_mut_dic_first['CG']/first_g.count)
                                
            if (current_type=="AC"):
                          prob_w=prob_w + (count_mut_dic_first['AC']/first_g.count) 
                          
            if (current_type=="TC"):
                          prob_w=prob_w + (count_mut_dic_first['TC']/first_g.count)  
                          
            if (current_type=="GC"):
                          prob_w=prob_w + (count_mut_dic_first['GC']/first_g.count) 
                          
            if (current_type=="AT"):
                          prob_w=prob_w + (count_mut_dic_first['AT']/first_g.count) 
                              
# =============================================================================
#  Finding the sum pf probabilities for second gene (p(cw) in (first,sec) order)  
# =============================================================================
           
            #count_mut_dic_sec=count_mutationType(sec_g) 
            #for t in range (len(sec_g_neo_mut)):                                    
            count_mut_dic_sec=count_mutationNeo(sec_g_neo_mut)
              
            prob_cw=0                # probability for all mutations for first gene in (first,sec) order
            for j in range(cw_len):
                current_type=cw_sub[j]
                if (current_type=="CT"):
                    prob_cw=prob_cw + (count_mut_dic_sec['CT']/sec_g.count)
                    
                if (current_type=="TA"):
                        prob_cw=prob_cw + (count_mut_dic_sec['TA']/sec_g.count)
                        
                if (current_type=="AG"):
                              prob_cw=prob_cw + (count_mut_dic_sec['AG']/sec_g.count)
                              
                if (current_type=="CA"):
                              prob_cw=prob_cw + (count_mut_dic_sec['CA']/sec_g.count)
                                    
                if (current_type=="GT"):
                              prob_cw=prob_cw + (count_mut_dic_sec['GT']/sec_g.count) 
                              
                if (current_type=="GA"):
                              prob_cw=prob_cw + (count_mut_dic_sec['GA']/sec_g.count)
                        
                if (current_type=="TG"):
                              prob_cw=prob_cw + (count_mut_dic_sec['TG']/sec_g.count)
                              
                if (current_type=="CG"):
                              prob_cw=prob_cw + (count_mut_dic_sec['CG']/sec_g.count)
                                    
                if (current_type=="AC"):
                              prob_cw=prob_cw + (count_mut_dic_sec['AC']/sec_g.count) 
                              
                if (current_type=="TC"):
                              prob_cw=prob_cw + (count_mut_dic_sec['TC']/sec_g.count)              
                              
                if (current_type=="GC"):
                              prob_cw=prob_cw + (count_mut_dic_sec['GC']/sec_g.count) 
                              
                if (current_type=="AT"):
                              prob_cw=prob_cw + (count_mut_dic_sec['AT']/sec_g.count)     
                              
       #['CT','TA','AG','CA','GT','GA','TG','CG','AC','TC','GC','AT']                        

# =============================================================================
# # Overall  similarity based on type between first and second gene 
# =============================================================================
            sigma_type= prob_w +prob_cw - 1         # sigma= p(w) +p(cw) - 1 
            
     else:
         sigma_type=1
       
     return sigma_type  



# ****************************************************************************************************
# ****************************************************************************************************
# ******************************************Main******************************************************
# ****************************************************************************************************
# ****************************************************************************************************



df1 = pd.read_spss('C:/GeneSimNet/OV_TP_MutationData_SNP_Full.sav')
arr=df1.to_numpy()    
row,col=arr.shape    # 5548*98
data=[]


# create a list of unique genes

for sample in range (row):
    if (sample==0):
        
        index=0
        count=1
        
        # read data related to the first gene
        current_sample_name=arr[sample][0]
        
        temp=Gene()
        temp.setName(current_sample_name)
        temp.setNeo(arr[sample][5])
        temp.setType(arr[sample][3])
        temp.setIndex(index)
        
        #if this is the only gene add it to list of genes
        if (sample == len(data)-1):
            temp.setCount(count)
            data.append(temp)
            
    else: # sample is not the first element  
     
        # if it is a repeated gene  
        if (current_sample_name == arr[sample][0]):
        
            temp.setNeo(arr[sample][5])
            temp.setType(arr[sample][3])    # NCBI_Build
            count=count+1
         
        
            #add last gene to list of genes
            if (sample == len(data)-1):
                data.append(temp)
      
        else:   # new gene
     
            # add previous graph to list of graphs
            if (temp.loaded):
                temp.setCount(count)
                data.append(temp)
                
            count=0
            #load new gene starting with the first connection
            temp = Gene() 
        
            current_sample_name = arr[sample][0]
            index = index + 1
            count=count+1
              
            temp.setNeo(arr[sample][5])
            temp.setType(arr[sample][3])
            temp.setName(arr[sample][0])
            temp.setIndex(index)
 
# =============================================================================
#  #Computing the similarity for type  and neoclotides

# =============================================================================
all_types=['Missense_Mutation','Nonsense_Mutation','RNA','Silent','Splice_Site','Nonstop_Mutation'] 
all_neos=['CT','TA','AG','CA','GT','GA','TG','CG','AC','TC','GC','AT']  
  



sim_matrix_type=np.zeros((len (data), len (data)))

for i in range (len(data)):
    for k in range (len(data)):
     if i==k:
        sim_matrix_type[i,k]=0    # if probability is 1 , there is no commonality ; if it is 0 they are very similar
     elif i!=k:
         first=data[i]       # first object which is the first gene
         sec=data[k]       # second object which is the second gene

    
#  Part 1-  calculating similarity based on type
         sigma_type_first_sec= similarity_type(first,sec)      # sigma (first, sec)= p(w)
         sigma_type_sec_first= similarity_type(sec,first)      # sigma (first, sec)= p(w)
 
         # find the max          
         temp = max(sigma_type_first_sec,sigma_type_sec_first)
         sim_matrix_type[i,k]=temp
     
 #  Part 2-  calculating similarity based on neocleotides  
         sigma_neo_first_sec= similarity_neo(first,sec)      # sigma (first, sec)= p(w)
         sigma_neo_sec_first= similarity_neo(sec,first)      # sigma (first, sec)= p(w)
         temp= max(sigma_neo_first_sec,sigma_neo_sec_first)
         sim_matrix_type[i,k]=(sim_matrix_type[i,k]+temp)/2.0       # add the similarity of neo to the similarity of type
         
         if (sim_matrix_type[i,k])<0 :
             print('i: ', i,'k:', k, 'value: ', sim_matrix_type[i,k] )

##############################################

# Get the name of all the genes from data  
gene_names_list=[]

for i in range (len(data)):
    gene_names_list.append(data[i].name)

gene_names = np.asarray(gene_names_list)