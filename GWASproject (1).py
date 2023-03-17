#!/usr/bin/env python
# coding: utf-8

# In[4]:


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import urllib.request
import gwaslab as gl
import pandas as pd
import argparse
import sys
import os

def process_paths(input_path,output_path):
    #Code to process each of the paths
    print(""" 
   ██████  ██     ██  █████  ███████     ██████  ██████   ██████       ██ ███████  ██████ ████████ 
  ██       ██     ██ ██   ██ ██          ██   ██ ██   ██ ██    ██      ██ ██      ██         ██    
  ██   ███ ██  █  ██ ███████ ███████     ██████  ██████  ██    ██      ██ █████   ██         ██    
  ██    ██ ██ ███ ██ ██   ██      ██     ██      ██   ██ ██    ██ ██   ██ ██      ██         ██    
   ██████   ███ ███  ██   ██ ███████     ██      ██   ██  ██████   █████  ███████  ██████    ██   
    """)
    print('Processing the paths...')
    return input_path, output_path

def process_data(input_path,output_path,cut,sig_level,skip,highlight,pinpoint,pinpoint_color, anno,highlight_color,build,xymt):
    
   # def download_ref_file(build):
    #    if build not in ["19","38"]:
    #        raise ValueError("Invalid build parameter")     
            
    #filename = f"refseq_hg{build}_gtf"

    i = 0
    for input_path in input_path:
        df=pd.read_csv(input_path,sep="\t",nrows=1000000)
        print('File data', input_path)
        print(df.head())
        
        #gl.check_available_ref()
        #gl.download_ref(filename)
        
        mysumstats = gl.Sumstats(
            df,
            fmt="plink",
            #For some reason, gwaslab understands that RSID = SNPID,
            #but it follows the same logic that SNPID have the format
            #chr:pos:ea:nea and RSID have the format 'rs1234'
            snpid="variant_id",
            #rsid="variant_id",
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
            build=build
        )
        mysumstats.infer_build()
       
        mysumstats.basic_check()
        
        mysumstats.get_lead(
           sig_level=sig_level,
           xymt=xymt,
           anno=anno,
           build=build,
           source="ensembl",
           verbose=True)
    
        file_name = output_path + 'Manhattan' + str(i) +'.pdf'
        #file_name2 = output_path + 'Regional' + str(i) +'.pdf'
        #file_name3 = output_path + 'Regional_other' + str(i) +'.pdf'  
        
        mysumstats.plot_mqq(save=file_name,
            saveargs={"dpi":400,"facecolor":"white"},
            cut=cut,
            sig_level=sig_level,
            anno=anno,
            skip=skip,
            highlight=highlight,
            pinpoint=pinpoint,
            pinpoint_color=pinpoint_color,
            highlight_color=highlight_color              
        )
             
    #mysumstats.get_lead(anno=anno)
   # mysumstats.plot_mqq(save= file_name2,region=(7,156538803,157538803),vcf_path=gl.get_path(filename))
        #mysumstats.plot_mqq(save= file_name3, mode="r", region=(7,156538803,157538803),region_grid=True, #gtf_path="ensembl")
        
    i +=1

#creating the parser
parser = argparse.ArgumentParser(description='Process multiple input paths')

#Adding the input path argument
parser.add_argument('--input_path', type=str, nargs='+', help='input paths')
parser.add_argument('--output_path', type=str, help='output path')
parser.add_argument('--sig_level', type=float,default=5e-8,help='sig_level value for plot_mqq')
parser.add_argument('--cut', type=int, default=0, help='cut value for plot_mqq')
parser.add_argument('--skip', type=int, default=0, help='skip value for plot_mqq')
parser.add_argument('--highlight',default=[],nargs="*", help='highlight value for plot_mqq')
parser.add_argument('--pinpoint', type=list, default=[], help='pinpoint value for plot_mqq')
parser.add_argument('--pinpoint_color', type=str, default="red", help='pinpoint_color value for plot_mqq')
parser.add_argument('--anno',type=bool,default="GENENAME",help='The variants to annotate will be selected automatically using a sliding window with windowsize=500kb')
parser.add_argument('--highlight_color', type=str, default='#CB132D', help='highlight_color value for plot_mqq')
parser.add_argument('--build',default=None,help='')
parser.add_argument('--xymt',default=True,help='')


#parsing the arguments
args = parser.parse_args()

#processing the paths
input_path, output_path = process_paths(args.input_path, args.output_path)

#processing the data
process_data(args.input_path, args.output_path, args.cut,args.sig_level, args.skip,args.highlight, args.pinpoint, args.pinpoint_color, args.anno,args.highlight_color,args.build,args.xymt)

# In[ ]:




