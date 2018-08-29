import os
import numpy as np
import bz2
import argparse

# uses similar distance measure as pyclone

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


def comb_list(seq, k):
    "returns a list of all k-combinations of the elements of sequence seq"
    n=len(seq)
    if not 0<=k<=n:
        raise Exception,"0<=k<=len(seq) is not true"
    v=[]   #list of combinations

    def f(x,y,a):
        if x==k:
            #we have taken enough elements, reject all remaining elements
            v.append(a)
            return
        if y==n-k:
            #we have rejected enough elements, take all remaining elements
            a.extend(seq[x+y:])
            v.append(a)
            return
        if (x<k):
            #take element seq[x+y]
            h=a+[seq[x+y]]
            f(x+1,y,h)
        if (y<n-k):
            #don't take element seq[x+y]
            f(x,y+1,a)            
    f(0,0,[])
    #print len(v)
    return v



def find_sample_passage_combinations(working_dir,trace_dir):
    
    stages_list=[]
    root_trace_path=working_dir+'/'+trace_dir
    for file in os.listdir(root_trace_path):
        if file.find('cellular_frequencies.tsv.bz2')!=-1:
            stage=file.split('.')[0]
            stages_list.append(stage)
            
    stage_cmobinations=comb_list(stages_list, 2)
    
    return stage_cmobinations
    
    
    
def prevalence_reader(infile_path,burn_num,statistic='mean'):
    """ reads the output of MCMC runs in pyclone. 
    Returns a dictionary with SNP id as key and average prevalence as value. """
    # output looks like this {'a': 1.3333333333333333, 'c': 5.666666666666667, 'b': 5.0, 'd': 1.3333333333333333}
    res_dic={}
    std_dict={}
    snp_dist_dict={}
    
    infile=bz2.BZ2File(infile_path)
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
        if statistic == 'mean':
            res_dic[snp_ids[i]]=np.mean(all_data_arr[:,i])
        if statistic == 'median':
            res_dic[snp_ids[i]]=np.median(all_data_arr[:,i])
        # getting std
        std_dict[snp_ids[i]]=np.std(all_data_arr[:,i])
    
    return snp_dist_dict



def cluster_file_reader(file_path):
    # reads the cluster IDs from pyclone and store them in a dict
    res_dict={}
    infile=open(file_path)    
    
    l=infile.readlines()
    for line in l[1:]:
        tmp=line.strip().split('\t')
        if tmp[1] in res_dict.keys():
            res_dict[tmp[1]].add(tmp[0])
        else:
            res_dict[tmp[1]]=set([tmp[0]])
    infile.close()
    return res_dict


    
def output_writer(exp_no, test_statistics_dict,output_folder, start_point,end_point,folder_name):
    
    
    outfile_statistics=open(output_folder+'/'+exp_no+'_statistic__'+folder_name+'__'+start_point+'_'+end_point+'.csv','w')
    for k,v in test_statistics_dict.iteritems():
        outfile_statistics.write(str(k)+'\t'+str(v)+'\n')
    
    outfile_statistics.close()
    
    
def cluster_avg_dist_calculator(pyclone_cluster_dict, snps_dist_dict,burnin):
    # calculates a single average distribution for each cluster
    clust_dist_dict={}
    for cluster_id,members_set in pyclone_cluster_dict.iteritems():
        avg_arr=np.zeros(burnin)
        for snp_id in members_set:
            avg_arr+=snps_dist_dict[snp_id]
        avg_arr=avg_arr/float(len(members_set))
        clust_dist_dict[cluster_id]=avg_arr
    return clust_dist_dict




def credible_interval_writer(clusters_credibility_dict,output_folder, start_stage,end_stage,start_cluster_credible_inerval_dict,end_cluster_credible_inerval_dict, start_clust_prev_dict, end_cluster_prev_dict):
    outfile=open(output_folder+'/'+'credible_interval__'+start_stage+'__'+end_stage+'.csv','w')
    outfile.write('cluster id'+'\t'+'x1'+'\t'+'x2'+'\t'+'y1'+'\t'+'y2'+'\t'+
                  'start prev'+'\t'+'end prev'+'\t'+'cred. test result'+'\n')
    for cluster_id in clusters_credibility_dict.keys():
        x1=start_cluster_credible_inerval_dict[cluster_id][0]
        x2=start_cluster_credible_inerval_dict[cluster_id][1]
        y1=end_cluster_credible_inerval_dict[cluster_id][0]
        y2=end_cluster_credible_inerval_dict[cluster_id][1]
        start_prev=np.mean(start_clust_prev_dict[cluster_id])
        end_prev=np.mean(end_cluster_prev_dict[cluster_id])
        test_res= clusters_credibility_dict[cluster_id]
        
        outfile.write(str(cluster_id)+'\t'+str(x1)+'\t'+str(x2)+'\t'+str(y1)+'\t'+str(y2)+'\t'+
                      str(start_prev)+'\t'+str(end_prev)+'\t'+str(test_res)+'\n')
    outfile.close()
    

#def run(root_path, exp_no, start_point,end_point, output_folder , folder_name):  omentum_site_1.cellular_frequencies.tsv.bz2
def run(pyclone_working_dir,trace_dir, pyclone_cluster_path, start_stage, end_stage, output_folder,burnin,CI_len):    
    
    start_cluster_credible_inerval_dict={}
    end_cluster_credible_inerval_dict={}
    clusters_credibility_dict={}

    start_pyclone_posterior_path=pyclone_working_dir+'/'+trace_dir+'/'+start_stage+'.cellular_frequencies.tsv.bz2'
    end_pyclone_posterior_path=pyclone_working_dir+'/'+trace_dir+'/'+end_stage+'.cellular_frequencies.tsv.bz2'
    
    pyclone_cluster_dict=cluster_file_reader(pyclone_cluster_path)
    
    # reading pyclone prevalence data for each snp at the starting point
    start_snps_dist_dict = prevalence_reader(start_pyclone_posterior_path,burnin,'mean')
    end_snps_dist_dict = prevalence_reader(end_pyclone_posterior_path,burnin,'mean')
    
    # calculating single distribution of clusters at start and end points
    start_clusters_dist_dict = cluster_avg_dist_calculator(pyclone_cluster_dict, start_snps_dist_dict,burnin)
    end_clusters_dist_dict = cluster_avg_dist_calculator(pyclone_cluster_dict, end_snps_dist_dict,burnin)
    
    a=float(100-CI_len)/2
    start=a
    end=100-a
    # calculating p-values for each cluster
    for cluster_id, snp_set in pyclone_cluster_dict.iteritems():
        # concatenating for start
        tmp=np.array([])
        for snp in snp_set:
            tmp=np.append(tmp,start_snps_dist_dict[snp])
        start_quantiles=np.percentile(tmp,[start,end])
        start_cluster_credible_inerval_dict[cluster_id]=start_quantiles
        
        # concatenating for end
        tmp=np.array([])
        for snp in snp_set:
            tmp=np.append(tmp,end_snps_dist_dict[snp])
        end_quantiles=np.percentile(tmp,[start,end])
        end_cluster_credible_inerval_dict[cluster_id]=end_quantiles
        
    # measuring the credibility of the passage
    for cluster_id in pyclone_cluster_dict.keys():
        x1= start_cluster_credible_inerval_dict[cluster_id][0]
        x2=start_cluster_credible_inerval_dict[cluster_id][1]
        y1=end_cluster_credible_inerval_dict[cluster_id][0]
        y2=end_cluster_credible_inerval_dict[cluster_id][1]
        #### calculating the difference
        prev_diff = abs(np.mean(start_clusters_dist_dict[cluster_id]) - np.mean(end_clusters_dist_dict[cluster_id]))

        threshd=0
        if (abs(x1-y2+threshd) >= abs(abs(x1-x2) + abs(y1-y2)) or abs(y1-x2+threshd) >= abs(abs(x1-x2) + abs(y1-y2))) and prev_diff >= 0.04 :
            clusters_credibility_dict[cluster_id]=0
        else:
            clusters_credibility_dict[cluster_id]=1
            
    credible_interval_writer(clusters_credibility_dict,output_folder, start_stage,end_stage,start_cluster_credible_inerval_dict,end_cluster_credible_inerval_dict, start_clusters_dist_dict, end_clusters_dist_dict)
    


if __name__=='__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="PyClone config file path")
    parser.add_argument("-f", help="path to cluster file.")
    parser.add_argument("-o", help="path to output folder.")
    parser.add_argument('-b',help="burnin")
    parser.add_argument("-i",help="length of the credible interval")
    
    
    args = parser.parse_args()
    
    confige_file_path=args.c
    pyclone_cluster_path=args.f
    output_folder=args.o
    burnin=int(args.b)
    CI_len=float(args.i)
    
    config_dict=read_config_file(confige_file_path)
    stage_combinations=find_sample_passage_combinations(config_dict['working_dir'], config_dict['trace_dir'])
    for stage_comb in stage_combinations:
        run(config_dict['working_dir'],config_dict['trace_dir'], pyclone_cluster_path, stage_comb[0], stage_comb[1], output_folder,burnin,CI_len)
    
    

    
    
    
    
    
    

    









