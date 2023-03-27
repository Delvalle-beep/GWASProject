#!/usr/bin/env python
# coding: utf-8

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import gwaslab as gl
import pandas as pd
import argparse
import sys
import os


def process_data(
    input_path,
    output_path,
    cut,
    sig_level,
    skip,
    highlight,
    pinpoint,
    pinpoint_color,
    anno,
    highlight_color,
    build,
    xymt,
    chr_filter
):
    
    def download_ref_file(build):
        if build not in ["19","38"]:
            raise ValueError("Invalid build parameter")     
            
    filename = f"ensembl_hg{build}_gtf"

    i = 0
    for input_path in input_path:
        df=pd.read_csv(input_path,sep="\t")
        print('File data', input_path)
        if 'snpid' or 'SNPID' not in df.columns:
            df['snpid'] = df.apply(lambda row: f"{row['chromosome']}:
                                   {row['base_pair_location']}:
                                   {row['effect_allele']}:
                                   {row['other_allele']}", 
                                   axis=1)
     
        gl.download_ref(filename)
#---------------------------------------------------------------------------------------------# 
        ##ISSUE##
        gtf_df = pd.read_csv(gl.get_path(filename),
                             sep="\t",
                             dtype={"gene_id": str,
                                    "attribute": str},
                             low_memory=False,
                             usecols=[
                                 0,
                                 1,
                                 2,
                                 3,
                                 4,
                                 5,
                                 6,
                                 7,
                                 8], 
                             names=["seqname",
                                    "source",
                                    "feature",
                                    "start",
                                    "end",
                                    "score",
                                    "strand",
                                    "frame",
                                    "attribute"]
                            )
        gtf_df.rename(columns={'seqname': 'chromosome'}, inplace=True)
        gtf_df['gene_name'] = gtf_df['attribute'].str.extract('gene_name "([^"]+)"', expand=False)
        gtf_df['chromosome'] = gtf_df['chromosome'].replace({'X': '23', 'Y': '24', 'MT': '25'})
        gtf_df['chromosome']=gtf_df['chromosome'].astype(int)
        
        merge_df= pd.merge(df,gtf_df,how='left',on='chromosome')

        #I can't merge, the error says "Unable to allocate 2.71 TiB for an array with shape (373042411878,) and data type int64", and I can't fix it even doing the data conversion
 #---------------------------------------------------------------------------------------------#       
        mysumstats = gl.Sumstats(
            merged_df,
            fmt="plink",
            snpid="snpid",
            rsid="variant_id",
            chrom="chromosome",
            pos="base_pair_location",
            ea="effect_allele",
            nea="other_allele",
            beta="beta",
            se="standard_error",
            p="p_value",
            build=build
        )
        
        mysumstats.infer_build()
        
        mysumstats.fix_id(
            fixchrpos=True
        )
        mysumstats.fix_chr(
            x=(23,"X"),
            y=(24,"Y"),
            mt=(25,"MT")
        )
        mysumstats.basic_check()
        
        mysumstats.assign_rsid(
            ref_rsid_tsv = gl.get_path(filename),
            chr_dict = gl.get_number_to_NC(build=build)
        )
        
        mysumstats.get_lead(
           sig_level=sig_level,
           xymt=xymt,
           anno=anno,
           source="ensembl",
           verbose=True
        )
    
        file_name = output_path + 'Manhattan' + str(i) +'.png'
        #file_name2 = output_path + 'Regional' + str(i) +'.pdf'
        #file_name3 = output_path + 'Regional_other' + str(i) +'.pdf'  
        
        if args.chr_filter:
            mysumstats.filter_value('CHR>1 & CHR<5').plot_mqq(
                save=file_name,
                saveargs={"dpi":100,"facecolor":"white"},
                cut=cut,
                sig_level=sig_level,
                anno=anno,
                skip=skip,
                highlight=highlight,
                pinpoint=pinpoint,
                pinpoint_color=pinpoint_color,
                highlight_color=highlight_color              
            )
        else:
            mysumstats.plot_mqq(
                save=file_name,
                saveargs={"dpi":100,"facecolor":"white"},
                cut=cut,
                sig_level=sig_level,
                anno=anno,
                skip=skip,
                highlight=highlight,
                pinpoint=pinpoint,
                pinpoint_color=pinpoint_color,
                highlight_color=highlight_color              
            )
             
      # mysumstats.get_lead(anno=anno)
      # mysumstats.plot_mqq(save= file_name2,region=(7,156538803,157538803),vcf_path=gl.get_path(filename))
      # mysumstats.plot_mqq(save= file_name3, mode="r", region=(7,156538803,157538803),region_grid=True, #gtf_path="ensembl")
        
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
parser.add_argument('--anno',type=bool,default=True,help='The variants to annotate will be selected automatically using a sliding window with windowsize=500kb')
parser.add_argument('--highlight_color', type=str, default='#CB132D', help='highlight_color value for plot_mqq')
parser.add_argument('--build',help='')
parser.add_argument('--xymt',default=True,help='')
parser.add_argument('--chr_filter',help='')

#parsing the arguments
args = parser.parse_args()

#processing the data
process_data(args.input_path,
             args.output_path,
             args.cut,
             args.sig_level, 
             args.skip,
             args.highlight, 
             args.pinpoint, 
             args.pinpoint_color, 
             args.anno,
             args.highlight_color,
             args.build,
             args.xymt,
             args.chr_filter)




