# updates yaml and posterior density files from pyclone output by removing the SNPs that break the credible interval

import os
import sys
import math
import numpy as np
import bz2
import pandas as pd
import copy
from shutil import copyfileobj
import yaml
import argparse


def read_config_file(confige_file_path):
    infile=open(confige_file_path)
    config_dict={}
    yaml_files_list=[]
    for line in infile:
        if line.find('num_iters')!=-1:
            tmp=line.strip().split(':')
            iter_num=int(tmp[1].strip())
            config_dict['num_iters']=iter_num
            
        if line.find('working_dir')!=-1:
            tmp=line.strip().split(':')
            working_dir_path=tmp[1].strip()
            config_dict['working_dir']=working_dir_path
            
        if line.find('trace_dir')!=-1:
            tmp=line.strip().split(':')
            trace_path=tmp[1].strip()
            config_dict['trace_dir']=trace_path
            
        if line.find('mutations_file')!=-1:
            tmp=line.strip().split(':')
            yaml_path=tmp[1].strip()
            yaml_files_list.append(yaml_path)
            config_dict['yaml_files_list']=yaml_files_list
                  
    return config_dict


def compress_file(infile_path,outfile_path):
    
    with open(infile_path, 'rb') as input:
        with bz2.BZ2File(outfile_path, 'wb', compresslevel=9) as output:
            copyfileobj(input, output)
            
def update_posterior_labels(infile_path, remove_list):
        
    infile=bz2.BZ2File(infile_path,'r')
    #infile=open(infile_path)
    snp_ids_list=infile.readline().strip().split('\t')
    
    my_data = np.genfromtxt(infile_path, skip_header=1, delimiter='\t')#,dtype=None)#, names=True)
    
    # making the list of indices of the SNPs to be removed from the posterior
    removed_snp_indx_list=[]
    for snp_id in remove_list:
        removed_snp_indx_list.append(snp_ids_list.index(snp_id))
    
    # deleting them from the posterior
    my_data = np.delete(my_data, removed_snp_indx_list,1)
    
    # updating the list of snp ids
    new_header_list=copy.deepcopy(snp_ids_list)
    for snp_id in remove_list:
        new_header_list.remove(snp_id)
    
    # Rename the previous files
    os.rename(infile_path, infile_path+'.original')
    
    df = pd.DataFrame(my_data,  columns=new_header_list)
    df.to_csv(infile_path[0:-4], index=False, header=True, sep='\t',float_format='%.0f')
    compress_file(infile_path[0:-4],infile_path)
    
    # removing the original file after the compression
    os.remove(infile_path[0:-4])
    
    

def update_posterior(infile_path, remove_list):
    
    if infile_path.find('bz2')!=-1:    
        infile=bz2.BZ2File(infile_path,'r')
    else:
        infile=open(infile_path)
    #infile=open(infile_path)
    snp_ids_list=infile.readline().strip().split('\t')
    
    my_data = np.genfromtxt(infile_path, skip_header=1, delimiter='\t')#,dtype=None)#, names=True)
    
    # making the list of indices of the SNPs to be removed from the posterior
    removed_snp_indx_list=[]
    for snp_id in remove_list:
        removed_snp_indx_list.append(snp_ids_list.index(snp_id))
    
    # deleting them from the posterior
    my_data = np.delete(my_data, removed_snp_indx_list,1)
    
    # updating the list of snp ids
    new_header_list=copy.deepcopy(snp_ids_list)
    for snp_id in remove_list:
        new_header_list.remove(snp_id)
    
    # Rename the previous files
    os.rename(infile_path, infile_path+'.original')

    
    df = pd.DataFrame(my_data,  columns=new_header_list)
    df.to_csv(infile_path[0:-4], index=False, header=True, sep='\t')
    compress_file(infile_path[0:-4],infile_path)
    
    # removing the original file after the compression
    os.remove(infile_path[0:-4])

        

def update_yaml(yaml_file_path,remove_list):
    
    removable_elements=[]
    with open(yaml_file_path,'r') as f:
        doc=yaml.load(f)
        for record_dict in doc['mutations']:
            if record_dict['id'] in remove_list:
                removable_elements.append(record_dict)
                
        for rec in removable_elements:
            doc['mutations'].remove(rec)
        
        # renaming the original file
    os.rename(yaml_file_path, yaml_file_path+'.original')
    stream=file(yaml_file_path,'w')
    yaml.dump(doc,stream)     
      
       
            
def update_all_files(root_pyclone_path,trace_dir,yaml_files_list,cluster_file_path,config_file_path, removable_snps_list,burnin):
    
    # updating lables.tsv
    root_trace_path=root_pyclone_path+'/'+trace_dir                    
    for file in os.listdir(root_trace_path):
        if file.find('labels')!=-1:
            labels_file_path=root_trace_path+'/'+file
            update_posterior_labels(labels_file_path,removable_snps_list)
            print 'updated'
                      
  
    # updating the yaml files 
      
    for yaml_file_path in yaml_files_list:
        update_yaml(root_pyclone_path+'/'+yaml_file_path,removable_snps_list)
        print 'update '+yaml_file_path
                     
                     
    # updating the posterior distributions
    for file in os.listdir(root_trace_path):
        if file.find('cellular_frequencies.tsv')!=-1:
            trace_file_path=root_trace_path+'/'+file
            update_posterior(trace_file_path,removable_snps_list)
            print 'update '+trace_file_path

                   
    # do the clustering  os.system('PyClone cluster -h')
    command1='PyClone cluster '+config_file_path+' '+cluster_file_path+' --burnin '+burnin
    os.system(command1)
     
    #updating the similarity matrix
    command2 = 'PyClone plot_similarity_matrix '+config_file_path+' '+root_pyclone_path+'/similarity_matrix.pdf'+' --burnin '+burnin
    os.system(command2)
            
                    
                    
def removable_SNPs_list_reader(removable_snps_path):
    res_list=[]
    infile=open(removable_snps_path)
    for line in infile:
        tmp=line.strip().split(',')
        res_list=tmp
    infile.close()
    return res_list
                  
 
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="PyClone config file path")
    parser.add_argument("-s", help="List of SNPs to be removed.")
    parser.add_argument("-f", help="path to the new cluster file.")
    parser.add_argument('-b',help="burnin")
    
    args = parser.parse_args()
    
    confige_file_path=args.c
    removable_snps_path=args.s
    cluster_file_path=args.f
    burnin=args.b
    
    config_dict=read_config_file(confige_file_path)
    
    removable_snps_list = removable_SNPs_list_reader(removable_snps_path)
    
    update_all_files(config_dict['working_dir'], config_dict['trace_dir'], config_dict['yaml_files_list'],cluster_file_path,confige_file_path, removable_snps_list,burnin)
    


    
    
    
    
    
    
    


    
    