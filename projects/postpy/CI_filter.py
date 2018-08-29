import os
import sys
import numpy as np
import bz2
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
    
def prevalence_reader(infile_path,burn_num, statistic='mean'):
    """ reads the output of MCMC runs in pyclone. 
    Returns a dictionary with SNP id as key and average prevalence as value. """
    res_dic={}
    std_dict={}
    snp_dist_dict={}
    
    if infile_path.find('frequencies.tsv.bz2')!=-1:
        infile=bz2.BZ2File(infile_path)
    else:
        infile=open(infile_path)
        
    data=infile.readlines()
    snp_ids=data[0].strip().split('\t')
    
    all_data=[]
    # ignores the burnin lines
    for i in range(burn_num+1,len(data)):
        tmp=data[i].strip().split('\t')
        all_data.append(tmp)
    
    all_data_arr=np.array(all_data,dtype=float)
    # assining list of values to snp dictionary
    
    for i in range(len(snp_ids)):
        snp_dist_dict[snp_ids[i]]=all_data_arr[:,i]
    
    return snp_dist_dict


def credible_interval_checker(snp_dist_dict,threshold,CI_len):
    # gets snp distribution dict and returns the list of SNPs with credible intervals larger than the threshold 
    accepted_snps_set=set([])
    rejected_snps_set=set([])
    
    a=float(100-CI_len)/2
    
    start=a
    end=100-a
    
    for snp_id, dist in snp_dist_dict.iteritems():
        #CI=np.percentile(dist,[5,95])
        CI=np.percentile(dist,[start,end])
        if CI[1]-CI[0] >= threshold:
            rejected_snps_set.add(snp_id)
        else:
            accepted_snps_set.add(snp_id)
            
    return accepted_snps_set,rejected_snps_set



def get_snps_list(root_folder_path,trace_dir,burnin):
    # gets the list of snps
    trace_folder=root_folder_path+'/'+trace_dir
    for file in os.listdir(trace_folder):
        if file.find('frequencies.tsv')!=-1:
            snps_dist_dict=prevalence_reader(trace_folder+'/'+file,burnin)
            break
    return snps_dist_dict.keys()
    

def snp_filter(root_folder_path,threshold,threshold_rejection_no,trace_dir,burnin,CI_len):

    res_dict={}
    culling_snps_list=[]

    # filters the SNPs based on the threshold 
    
    trace_folder=root_folder_path+'/'+trace_dir
    
    # file_no for measuring the threshold
    file_no=0
    
    #initializing
    number_of_rejections_per_snp_dict={}
    snp_id_list=get_snps_list(root_folder_path,trace_dir,burnin)
    for snp_id in snp_id_list:
        number_of_rejections_per_snp_dict[snp_id]=0

    for file in os.listdir(trace_folder):
        if file.find('frequencies.tsv')!=-1:
  
            
            # file_no for measuring the threshold
            file_no+=1
            
            stage=file.split('.')[0]
            snps_dist_dict=prevalence_reader(trace_folder+'/'+file,burnin)
            accepted_snps_set,rejected_snps_set = credible_interval_checker(snps_dist_dict,threshold,CI_len)
            res_dict[stage]=rejected_snps_set
                
            for snp_id in res_dict[stage]:
                number_of_rejections_per_snp_dict[snp_id]+=1
        
        #####
    for snp, rej_no in number_of_rejections_per_snp_dict.iteritems():
        if float(rej_no)/file_no > threshold_rejection_no:
            culling_snps_list.append(snp)           
        culling_snps_percentage=float(len(culling_snps_list))/len(snps_dist_dict.keys())    
               
    return culling_snps_list,culling_snps_percentage


    
    
def result_writer(culling_snps_list,culling_snps_percentage,outfile_path,threshold_rejection_no,threshold):
    # in culling_snps_percentage_dict the key is the folder name and the value is the percentage of the rejected mutations
    outfile=open(outfile_path+'/CI_'+str(threshold)+'_node_rej_'+str(threshold_rejection_no)+'.csv','w')
    
    outfile.write('Threshold rej., maximum CI, precentage of removed SNVs \n')
    outfile.write(str(threshold_rejection_no)+','+str(threshold)+','+str(culling_snps_percentage))
    outfile.close()
    
    # writing the list of the snps that should be culled
    snps_list_filepath=outfile_path+'/snps_list_CI_'+str(threshold)+'_node_rej_'+str(threshold_rejection_no)+'.csv'
    snps_listfile=open(snps_list_filepath,'w')

    tmp_s=''
    for snp_id in culling_snps_list:
        tmp_s+=','+snp_id
    snps_listfile.write(tmp_s[1:]+'\n')
    snps_listfile.close()
    
    
def snps_list_writer(culling_snps_list,culling_snps_percentage,outfile_path,threshold_rejection_no,threshold):
    # in culling_snps_percentage_dict the key is the folder name and the value is the percentage of the rejected mutations
    
    # writing the list of the snps that should be culled
    snps_list_filepath=outfile_path+'/snps_list_CI_'+str(threshold)+'_node_rej_'+str(threshold_rejection_no)+'.csv'
    snps_listfile=open(snps_list_filepath,'w')

    tmp_s=''
    for snp_id in culling_snps_list:
        tmp_s+=','+snp_id
    snps_listfile.write(tmp_s[1:]+'\n')
    snps_listfile.close()


def stat_writer(outfile_path,threshold_rejection_no,res_list):
    outfile=open(outfile_path+'/summary_'+str(threshold_rejection_no)+'.csv','w')
    outfile.write('Threshold rej., maximum CI, fraction of removed SNVs \n')
    for r in res_list:
        outfile.write(str(threshold_rejection_no)+','+str(r[0])+','+str(r[1])+'\n')
    outfile.close()
    
    
if __name__=='__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("-r", help="Rejection threshold")
    parser.add_argument("-c", help="PyClone config file path")
    parser.add_argument("-o", help="outputs path")
    parser.add_argument('-b',help="burnin")
    parser.add_argument("-i",help="length of the credible interval")

    
    args = parser.parse_args()
    if args.r is None:
        threshold_rejection_no = 0.2
    else:
        threshold_rejection_no =float(args.r)
        
    #root_folder_path =args.f
    outfile_path=args.o
    config_dict=read_config_file(args.c)
    burnin=int(args.b)
    CI_len=float(args.i)
    
    res_list=[]
    for threshold in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        culling_snps_list,culling_snps_percentage=snp_filter(config_dict['working_dir'], threshold, threshold_rejection_no,config_dict['trace_dir'],burnin,CI_len)
        print str(threshold)+'  done!'
        snps_list_writer(culling_snps_list, culling_snps_percentage,outfile_path,threshold_rejection_no,threshold)
        ##
        res_list.append([threshold,culling_snps_percentage])
    
    stat_writer(outfile_path,threshold_rejection_no,res_list)
    

    








