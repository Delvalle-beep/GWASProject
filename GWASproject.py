#!/usr/bin/env python
# coding: utf-8

# In[3]:


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import gwaslab as gl
import pandas as pd
import argparse
import sys
import os

def process_paths(input_path,output_path):
    #Code to process each of the paths
    print('Processing the paths...')
    return input_path, output_path

def process_data(input_path,output_path):
    i = 0
    for input_path in input_path:
        df=pd.read_csv(input_path,sep="\t",nrows=100000)
        print('File data', input_path)
        print(df.head())
        
        #loading the file in gwaslab
        mysumstats = gl.Sumstats(
             input_path,
             fmt="plink",
             snpid="variant_id",
             rsid=None,
             chrom="chromosome",
             pos="base_pair_location",
             ea="effect_allele",
             nea="other_allele",
             ref=None,
             alt=None,
             eaf=None,
             neaf=None,
             n=None,
             beta="beta",
             se="standard_error",
             chisq=None,
             z=None,
             p="p_value",
             mlog10p=None,
             info=None,
             OR="odds_ratio",
             OR_95L=None,
             OR_95U=None,
             status=None,
             other=[],
             direction=None,
             verbose=True,
             build="99"
            )
        
        mysumstats.random_variants(n=100000,inplace=True)
        mysumstats.basic_check()
        mysumstats.infer_build()
        mysumstats.lookup_status()
        
        file_name = output_path + 'Manhatan' + str(i) +'.pdf'
        file_name2 = output_path + 'Regional' + str(i) +'.pdf'
        file_name3 = output_path + 'Regional_other' + str(i) +'.pdf'
        mysumstats.plot_mqq(save= file_name,saveargs={"dpi":400,"facecolor":"white"})
        mysumstats.get_lead()
        mysumstats.plot_mqq(save= file_name2,region=(7,156538803,157538803))
        mysumstats.plot_mqq(save= file_name3, mode="r", region=(7,156538803,157538803),region_grid=True, gtf_path="ensembl")
        i +=1
    
    
#creating the parser
parser = argparse.ArgumentParser(description='Process multiple input paths')

#Adding the input path argument
parser.add_argument('--input_path', type=str, nargs='+', help='input path')
parser.add_argument('--output_path', type=str, help='output path')

#parsing the argument
args = parser.parse_args()

#using the arguments
input_path, output_path = process_paths(args.input_path, args.output_path)
process_paths(args.input_path,args.output_path)
                    





# In[18]:





# In[ ]:




