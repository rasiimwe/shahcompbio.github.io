"""
Created on Wed Oct 16 11:31:56 2013

@author: jtaghiyar
"""
from __future__ import division
import numpy
from math import log
from scipy.stats import binom
from scipy.special import gammaln


class Features:
    def __init__(self, tumour_tuple=None, normal_tuple=None, reference_tuple=None, purity=70):
        self.tt = tumour_tuple
        self.nt = normal_tuple
        self.rt = reference_tuple
        self.name = "TCGA Benchmark 4 feature set with coverage info"
        self.version = "4.1.2"
 
        ## check for zero coverage or None tuple
        if self.tt is None or self.tt[5][0] == 0:
            self.tt = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, 1, 1, [1], None)

        if self.nt is None or self.nt[5][0] == 0: 
            self.nt = (None, [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1,1,1,1,1,[1,1]], [1]*6, 1, 1, 1, 1, 1, 1, 1, 1, [1], None)
        
        if self.rt is None:
            self.rt = (0, 0, 0, 0, 0)
         
        ## reference base index + 1 = index of the same base in the tumour/normal bam tuple
        self.b  = self.rt[0] + 1  
        
        ## coverage data        
        self.cd = (float(30), float(30), int(purity), float(0))
        
        ## to avoid division by zero
        self.ep = 1e-5 

        ## get the the count of nonref bases from both tumour/normal tuples. Get only for the variant that is common between both and has the maximum depth.
        self.__get_tumour_nonref_index()
        self.__get_tumour_nonref_count()
        self.__get_normal_nonref_count()

        ## calculate jointsnvmix values
        #self.__jointsnvmix()
        self.__multinom_sb()

        ## expectation of success for binomial test
#         self.p = 0.01
        
#=============================================================================
# features
#=============================================================================
        self.coverage_features = (
        #enet:("normal_depth_coverage", self.nt[5][0] / self.cd[0]),
        #enet:("tumour_depth_coverage", self.tt[5][0] / self.cd[1]),
        #enet:("normal_contamination", self.cd[2] / 100),
        #enet:("whole_genome", self.cd[3])
        )
        
        self.feature_set = (
        ## log scale
        #("tumour_ref_depth_log", log(self.tt[self.b][0] / self.tt[5][0] + self.ep)),
        #("normal_ref_depth_log", log(self.nt[self.b][0] / self.nt[5][0] + self.ep)),
        
        ("homopolymer_b", self.rt[2]),
        ("homopolymer_f", self.rt[1]),
        ("normal_depth", self.nt[5][0]),
        #enet:("normal_direction_ratio", self.nt[5][4] / self.nt[5][0]),
        ("normal_distance_ratio", self.nt[5][3] / self.nt[5][0]),
        ("normal_entropy", self.nt[10]),
        #enet:("normal_indels", (self.nt[9]+self.nt[11]) / self.nt[5][0]),
        ("normal_mapq_ratio", self.nt[5][2] / self.nt[5][0]),
        ("normal_quality_ratio", self.nt[5][1] / self.nt[5][0]),
        #enet:("normal_ref_depth", self.nt[self.b][0] / self.nt[5][0]),
        #enet:("normal_ref_direction", self.nt[self.b][4] / (self.nt[self.b][0] + self.ep)),
        #enet:("normal_ref_direction_total", self.nt[self.b][4] / self.nt[5][0]),
        ("normal_ref_quality", self.nt[self.b][1] / (self.nt[self.b][0] + self.ep)), #self.nt[5][0]),
        ("normal_tumour_depth", self.tt[5][0] / self.nt[5][0]),
        ("normal_tumour_direction", (self.tt[5][4] / self.tt[5][0]) / ((self.nt[5][4] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_entropy", self.nt[10] / (self.tt[10] + 1e-8)),
        ("normal_tumour_mapq", (self.tt[5][2] / self.tt[5][0]) / ((self.nt[5][2] / self.nt[5][0]) + self.ep)),
        #enet:("normal_tumour_quality", (self.tt[5][1] / self.tt[5][0]) / ((self.nt[5][1] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_ref_depth", ((self.tt[self.b][0] / self.tt[5][0]) + self.ep) / ((self.nt[self.b][0] / self.nt[5][0]) + self.ep)),
        ("normal_tumour_ref_direction", (self.tt[self.b][4] / (self.tt[self.b][0] + self.ep)) / ((self.nt[self.b][4] / (self.nt[self.b][0] + self.ep)) + self.ep)),
        #enet:("normal_variant_allele_frequency", self.nt_nonref_count / self.nt[5][0]),
        #enet:("normal_variant_depth_ratio", self.nt_nonref_count / (self.nt[5][0] + self.ep)),
        #enet:("normal_variant_direction_ratio", self.nt[self.nonref_index][4] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_distance", self.nt[self.nonref_index][3] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_mapq_mean", self.nt[self.nonref_index][2] / (self.nt_nonref_count + self.ep)),
        ("normal_variant_quality_ratio", self.nt[self.nonref_index][1] / (self.nt_nonref_count + self.ep)),
        #enet:("region_entropy", self.rt[4]),
        #enet:("region_gc_content", self.rt[3]),
        ("tumour_depth", self.tt[5][0]),
        #enet:("tumour_direction_ratio", self.tt[5][4] / self.tt[5][0]),
        ("tumour_distance_ratio", self.tt[5][3] / self.tt[5][0]),
        #enet:("tumour_entropy", self.tt[10]),
        #enet:("tumour_indels", (self.tt[9]+self.tt[11]) / self.tt[5][0]),
        ("tumour_mapq_ratio", self.tt[5][2] / self.tt[5][0]),
        ("tumour_quality_ratio", self.tt[5][1] / self.tt[5][0]),
        #enet:("tumour_ref_depth", self.tt[self.b][0] / self.tt[5][0]),
        #enet:("tumour_ref_direction", self.tt[self.b][4] / (self.tt[self.b][0] + self.ep)),
        #enet:("tumour_ref_direction_total", self.tt[self.b][4] / self.tt[5][0]),
        ("tumour_ref_quality", self.tt[self.b][1] / (self.tt[self.b][0] + self.ep)), #self.tt[5][0]),
        #enet:("tumour_variant_allele_frequency", self.tt_nonref_count / self.tt[5][0]),
        #enet:("tumour_variant_depth_ratio", self.tt_nonref_count / (self.tt[5][0] + self.ep)),
        #enet:("tumour_variant_direction_ratio", self.tt[self.nonref_index][4] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_distance", self.tt[self.nonref_index][3] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_mapq_mean", self.tt[self.nonref_index][2] / (self.tt_nonref_count + self.ep)),
        ("tumour_variant_quality_ratio", self.tt[self.nonref_index][1] / (self.tt_nonref_count + self.ep)),

        #enet:("normal_direction_ratio_mean", self.nt[self.nonref_index][4]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        ("normal_mapq_ratio_mean", self.nt[self.nonref_index][2]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep) ),
        #enet:("normal_quality_ratio_mean", self.nt[self.nonref_index][1]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep) ),
        #enet:("normal_distance_ratio_mean", self.nt[self.nonref_index][3]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        #("normal_nonref_depth_mean", self.nt[self.nonref_index][0]/(self.nt[self.nonref_index][0]+self.ep) / (self.nt[5][0]+self.ep)),
        
        #enet:("tumour_direction_ratio_mean", self.tt[self.nonref_index][4]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        ("tumour_mapq_ratio_mean", self.tt[self.nonref_index][2]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        ("tumour_quality_ratio_mean", self.tt[self.nonref_index][1]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        #enet:("tumour_distance_ratio_mean", self.tt[self.nonref_index][3]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
        #("tumour_nonref_depth_mean", self.tt[self.nonref_index][0]/(self.tt[self.nonref_index][0]+self.ep) / (self.tt[5][0]+self.ep)),
                
        #("normal_tumour_depth_mean", ((self.nt[self.nonref_index][0]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][0]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]) ),
        ("normal_tumour_direction_mean", ((self.nt[self.nonref_index][4]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][4]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0] +self.ep) ),
        ("normal_tumour_mapq_mean",((self.nt[self.nonref_index][2]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][2]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
        ("normal_tumour_quality_mean", ((self.nt[self.nonref_index][1]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][1]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
        ("normal_tumour_distance_mean", ((self.nt[self.nonref_index][3]/(self.nt[self.nonref_index][0]+self.ep))/self.nt[5][0]) / ((self.tt[self.nonref_index][3]/(self.tt[self.nonref_index][0]+self.ep))/self.tt[5][0]+self.ep) ),
        
#         ## jointsnvmix
#         ("tumour_normal_jointsnvmix_somatic", self.p_somatic),
#         ("tumour_normal_jointsnvmix_germline", self.p_germline),
#         ("tumour_normal_jointsnvmix_wildtype", self.p_wildtype),
#         ("tumour_normal_jointsnvmix_loh", self.p_loh),
#         ("tumour_normal_jointsnvmix_error", self.p_error),

#         ("normal_tumour_variant_direction_ratio", (self.tt[self.nonref_index][4]) / (self.nt[self.nonref_index][4] + 1)), 
#         ("normal_tumour_variant_mapq_ratio", (self.tt[self.nonref_index][2] ) / (self.nt[self.nonref_index][2] + 1)), 
#         ("normal_tumour_variant_quality_ratio", (self.tt[self.nonref_index][1] ) / (self.nt[self.nonref_index][1] + 1)), 
#         ("normal_tumour_variant_distance_ratio", (self.tt[self.nonref_index][3] ) / (self.nt[self.nonref_index][3] + 1)), 
#         ("normal_tumour_distance", (self.tt[5][3] / self.tt[5][0]) / ((self.nt[5][3] / self.nt[5][0]) + self.ep)),

        ("normal_distance", self.nt[5][3]),
        ("normal_variant_direction", self.nt[self.nonref_index][4]),
        ("normal_variant_distance", self.nt[self.nonref_index][3]),
        ("normal_variant_mapq", self.nt[self.nonref_index][2]),
        ("normal_variant_quality", self.nt[self.nonref_index][1]),
        ("tumour_distance", self.tt[5][3]),
        ("tumour_variant_direction", self.tt[self.nonref_index][4]),
        ("tumour_variant_distance", self.tt[self.nonref_index][3]),
        ("tumour_variant_mapq", self.tt[self.nonref_index][2]),
        ("tumour_variant_quality", self.tt[self.nonref_index][1]),

        ("probability_pos", self.pos),
        #("normal_indel_ratio",sum(self.nt[-2])/self.nt[5][0] ),
        #("tumour_indel_ratio",sum(self.tt[-2])/self.tt[5][0] ),
        
#         ("normal_direction", self.nt[5][4]),
#         ("tumour_direction", self.tt[5][4]),
#         ("normal_variant_direction", self.nt[5][4] - self.nt[self.rt[0] + 1][4]),
#         ("tumour_variant_direction", self.tt[5][4] - self.tt[self.rt[0] + 1][4]),
#         ("normal_tumour_variant_ratio", (self.tt[5][4] - self.tt[self.rt[0] + 1][4]) / (self.nt[5][4] - self.nt[self.rt[0] + 1][4] + self.ep)),
#         ("normal_tumour_variant_quality", (self.tt[5][1] - self.tt[self.rt[0] + 1][1]) / (self.nt[5][1] - self.nt[self.rt[0] + 1][1] + self.ep)),
#         ("normal_tumour_variant_distance", (self.tt[5][3] - self.tt[self.rt[0] + 1][3]) / (self.nt[5][3] - self.nt[self.rt[0] + 1][3] + self.ep)),
#         ("normal_quality", self.nt[5][1]),
#         ("tumour_quality", self.tt[5][1]),
#         ("normal_tumour_variant_mapq", (self.tt[5][2] - self.tt[self.rt[0] + 1][2]) / (self.nt[5][2] - self.nt[self.rt[0] + 1][2] + self.ep)),
#         ("normal_minor_allele", (self.rt[0] != self.tt[7] and self.nt[self.tt[6] + 1][0] > 0) or (self.rt[0] != self.tt[6] and self.nt[self.tt[7] + 1][0] > 0)),
#         ("normal_variant_depth", (self.nt[5][0] - self.nt[self.b][0])),
#         ("tumour_variant_depth", (self.tt[5][0] - self.tt[self.b][0])),
#         ("tumour_mapq", self.tt[5][2]),
#         ("normal_mapq", self.nt[5][2]),
#         ("normal_variant_mapq", (self.nt[5][2] - self.nt[self.rt[0] + 1][2])),
#         ("tumour_variant_mapq", (self.tt[5][2] - self.tt[self.rt[0] + 1][2])),
       
        ## binomial test
#         ("tumour_variant_depth_binom", binom.pmf(self.tt_nonref_count, self.tt[5][0], self.p)),
#         ("normal_tumour_variant_depth", (self.tt[5][0] - self.tt[self.rt[0] + 1][0]) / (self.nt[5][0] - self.nt[self.rt[0] + 1][0] + self.ep)),
        
        ##TODO: should be fixed, changed to test
#         ("normal_tumour_variant_depth_ratio", min(((self.tt[5][0] - self.tt[self.b][0]) / self.tt[5][0]) / (((self.nt[5][0] - self.nt[self.b][0]) / self.nt[5][0]) + self.ep), 5)),
#         ("normal_tumour_variant_depth_ratio", (self.tt_nonref_count / self.tt[5][0]) / (self.nt_nonref_count / self.nt[5][0] + self.ep)),
#         ("normal_tumour_variant_depth_log_ratio", log((self.tt_nonref_count / self.tt[5][0]) + self.ep) / log(self.nt_nonref_count / self.nt[5][0] + self.ep)),

        )
        
    def __get_tumour_nonref_index(self):
        nonrefbases = [x for x in range(4) if x != self.b - 1]
        max_nonrefbase_depth = 0
        nonref_index = nonrefbases[1]
        if nonref_index == self.b:
            nonref_index = nonrefbases[2]  
        
        for nb in nonrefbases:
            index = nb + 1
            tumour_nonrefbase_depth = self.tt[index][0]
            if tumour_nonrefbase_depth != 0 and tumour_nonrefbase_depth > max_nonrefbase_depth:
                max_nonrefbase_depth = self.tt[index][0]
                nonref_index = index
                
        self.nonref_index = nonref_index
    
    def __get_tumour_nonref_count(self):
        self.tt_nonref_count = self.tt[self.nonref_index][0]

    def __get_normal_nonref_count(self):
        self.nt_nonref_count = self.nt[self.nonref_index][0]
    
    ##TODO: eventually this function should be removed    
    def __isvalid(self, x):
        if numpy.isnan(x) or numpy.isinf(x):
            
            ##TODO: remove this line
            print "NaN"
            return False
        return True
    
    def __xentropy(self):
        tumour_base_qualities = (self.tt[1][1], self.tt[2][1], self.tt[3][1], self.tt[4][1], self.tt[5][1])
        normal_base_qualities = (self.nt[1][1], self.nt[2][1], self.nt[3][1], self.nt[4][1], self.nt[5][1])
        total_tbq = tumour_base_qualities[4]
        total_nbq = normal_base_qualities[4]
        ent = 0 # entropy
        
        if total_tbq == 0 or total_nbq == 0:
            return ent
            
        for i in xrange(4):
            tumour_base_probability = tumour_base_qualities[i] / total_tbq
            normal_base_probability = normal_base_qualities[i] / total_nbq            
            if tumour_base_probability != 0:
                if normal_base_probability == 0:
                    ent -= -7 * tumour_base_probability
                else:
                    ent -= log(normal_base_probability) * tumour_base_probability
        return ent

    def __log_factorial(self, x):
        return gammaln(numpy.array(x)+1)
    
    def __multinomial(self, xs, ps):
        n = sum(xs)
        xs, ps = numpy.array(xs), numpy.array(ps)
        result = self.__log_factorial(n) - sum(self.__log_factorial(xs)) + sum(xs * numpy.log(ps))
        return numpy.exp(result)

    def __multinom_sb(self):
	t_ref_pos = self.tt[self.b][5][0]
        t_ref_neg = self.tt[self.b][5][1]
        t_nref_pos =self.tt[self.nonref_index][5][0] 
        t_nref_neg =self.tt[self.nonref_index][5][1]
        n_ref_pos = self.nt[self.b][5][0]
        n_ref_neg = self.nt[self.b][5][1]
        n_nref_pos =self.nt[self.nonref_index][5][0]
        n_nref_neg =self.nt[self.nonref_index][5][1]

        counts = [n_ref_pos,n_ref_neg,n_nref_pos,n_nref_neg,t_ref_pos,t_ref_neg,t_nref_pos,t_nref_neg]
        #prior
        ep = 0.001
        prob = [[1/6,1/6,ep,ep,1/6,1/6,1/6,1/6],
                [1/6,1/6,ep,ep,ep,ep,1/6,1/6],
                [1/6,1/6,ep,ep,1/6,1/6,0.1,0.1],
                [1/8,1/8,1/8,1/8,1/8,1/8,1/8,1/8],
                [ep,ep,1/4,1/4,ep,ep,1/4,1/4],
                [1/3,1/3,ep,ep,ep,ep,1/3,ep],
                [1/3,1/3,ep,ep,ep,ep,ep,1/3],
                [1/5,1/5,ep,ep,1/5,1/5,1/5,ep],
                [1/5,1/5,ep,ep,1/5,1/5,ep,1/5],
                [1/4,1/4,ep,ep,1/4,1/4,ep,ep]]

        prob = [[val/sum(prob[i]) for val in classval]for i,classval in enumerate(prob)]
        
        #output
        result = [self.__multinomial(counts, ps) for ps in prob]
        if sum(result) == 0:
            result = [val/(sum(result)+1e-300) for val in result]
        else:    
            result = [val/sum(result) for val in result]

        self.pos = result[0]+result[1]+result[2]
    
    def __jointsnvmix(self):
        ## a_priori
        aa_aa = 1e3
        aa_ab = 1e1
        aa_bb = 1e1
        ab_aa = 1e1
        ab_ab = 1e2
        ab_bb = 1e1
        bb_aa = 1
        bb_ab = 1
        bb_bb = 1e2
        a_priori_list = numpy.array([aa_aa, aa_ab, aa_bb, ab_aa, ab_ab, ab_bb, bb_aa, bb_ab, bb_bb])
        a_priori_mat  = a_priori_list.reshape(3,3)
        
        ## normal allele frequency
        aa = 0.01
        ab = 0.50
        bb = 0.99
        n_prob = [aa, ab, bb]
        
        ## tumour allele frequency
        aa = 0.01
        ab = 0.30
        bb = 0.90
        t_prob = [aa, ab, bb]
        
        ## binomial pmf for three different probabilities
        t_binom = [binom.pmf(self.tt_nonref_count, self.tt[5][0], p) for p in t_prob]
        n_binom = [binom.pmf(self.nt_nonref_count, self.nt[5][0], p) for p in n_prob]
        
        ## transpose(t_binom) x n_binome = 3x3 binom_mat matrix
        binom_mat = [[i*j for i in n_binom] for j in t_binom]
        binom_mat = numpy.array(binom_mat)
        
        ## element by element product of a_priori matrix to the m_binom matrix
        p_mat = numpy.multiply(a_priori_mat, binom_mat)
        
        ## normalize 
        weight = sum(sum(p_mat))
        p_mat = p_mat / (weight + self.ep)
    
        ## p_SOMATIC = p_aa,ab + p_aa,bb
        self.p_somatic = p_mat[1][0] + p_mat[2][0]
        
        ## p_WILDTYPE = p_aa,aa
        self.p_wildtype = p_mat[0][0]
        
        ## p_GERMLINE = p_ab,ab + p_bb,bb
        self.p_germline = p_mat[1][1] + p_mat[2][2]
        
        ## p_LOH = p_ab,aa + p_ab,bb
        self.p_loh = p_mat[0][1] + p_mat[2][1]
        
        ## p_ERROR = p_bb,aa + p_bb,ab
        self.p_error = p_mat[0][2] + p_mat[1][2]
        
    def get_features(self):
        features = []
        for _, f in self.feature_set:
            if self.__isvalid(f):
                features.append(f)
            else:
                features.append(0)
       
        for _,f in self.coverage_features:
            if self.__isvalid(f):
                features.append(f)
            else:
                features.append(0)
        
        features.append(self.__xentropy())
        return features
       
    def get_feature_names(self):
        feature_names = []
        for n, _ in self.feature_set:
            feature_names.append(n)

        for n, _ in self.coverage_features:
            feature_names.append(n)
            
        feature_names.append("xentropy")
        
        return feature_names








    
    
        
