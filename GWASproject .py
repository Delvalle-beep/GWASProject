#!/usr/bin/env python
# coding: utf-8
from matplotlib.backends.backend_pdf import PdfPages
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt
import urllib.request
import gwaslab as gl
import pandas as pd
import argparse
import sys
import os

def process_data(input_paths,output_path,cut,sig_level,skip,highlight,
                 pinpoint,pinpoint_color,anno,highlight_color,build,
                 chr_filter,vcf_file):
    
    #Checks if the build was declared a valid value
    if build not in ["19","38"]:
        raise ValueError("Invalid build parameter")
        
    #variable to define the reference file to be downloaded  
    filename = f"1kg_eas_hg{build}"
    
    #The counter has been defined to count the number of files received
    #in input and in this way save the plot according the file runs
    for i, input_path in enumerate(input_paths):
        df=pd.read_csv(input_path,sep="\t")
        df = df.dropna()
        print('File data', input_path)
        if 'snpid' or 'SNPID' not in df.columns:
            df['snpid'] = df.apply(lambda row: f"{row['chromosome']}:{row['base_pair_location']}:{row['effect_allele']}:{row['other_allele']}",axis=1)
        print(df.head())

        #If the user wants to include the LD in the plot 
        if args.vcf_file=='True':
            # Define the path where the file should be saved
            file_path = os.path.expanduser(f'~/.gwaslab/1kgp3v5.hg{build}.vcf.gz')
            # Check if the file already exists in the specified path
            if os.path.exists(file_path):
                print(f'File {filename} already exists at path {file_path}')
            else:
            # If the file does not exist, download it using the gl.download_ref() function
                print(f'The file {filename} was not found on the computer. Downloading the file...')
                gl.download_ref(filename)
                
        #If the user wants to include the LD from its own file
        elif args.vcf_file is not None:
            file_path = args.vcf_file
            if os.path.exists(file_path):
                print(f'The file {file_path} was found')
            else:
                print(f'The file {file_path}  was not found')
                sys.exit(1)
        #If the user doesn't set anything in the vcf_file, continue the script normally 
        else:
             pass

        #Loading the data and transforming it into PLINK format
        mysumstats = gl.Sumstats(
            df,
            fmt='plink',
            beta='beta',
            p='p_value',
            build=build,
            snpid='snpid',
            rsid='variant_id',
            ea='effect_allele',
            chrom='chromosome',
            nea='other_allele',
            se='standard_error',
            pos='base_pair_location',
        )

        mysumstats.basic_check()

        file_name = output_path + 'Manhattan' + str(i) +'.png'
        
        #file_name3 = output_path + 'Regional_other' + str(i) +'.pdf'

        #To filter the chromosomes or range to be plotted
        if args.chr_filter:
            mysumstats.filter_value(chr_filter).plot_mqq(
                mode="m",
                cut=cut,
                anno=anno,
                skip=skip,
                save=file_name,
                pinpoint=pinpoint,
                sig_level=sig_level,
                highlight=highlight,
                pinpoint_color=pinpoint_color,
                highlight_color=highlight_color,
                saveargs={"dpi":80,"facecolor":"white"}
            )
        else:
            mysumstats.plot_mqq(
                mode="m",
                cut=cut,
                anno=anno,
                skip=skip,
                save=file_name,
                pinpoint=pinpoint,
                sig_level=sig_level,
                highlight=highlight,
                pinpoint_color=pinpoint_color,
                highlight_color=highlight_color,
                saveargs={"dpi":80,"facecolor":"white"}
            )

        x = mysumstats.get_lead()
        print(x)
        
        for i, row in x.iterrows():
            file_name2 = output_path + 'Regional_Plot' + str(i) +'.png'
            chrom = row['CHR']
            pos = row['POS']
            region_start = pos - 500000
            region_end = pos + 500000
            
            if args.vcf_file=='True':
                mysumstats.plot_mqq(
                    save=file_name2,
                    vcf_path=gl.get_path(filename),
                    region=(chrom,region_start,region_end),
                    saveargs={"dpi":80,"facecolor":"white"}
                )
            elif args.vcf_file is not None:
                mysumstats.plot_mqq(
                    save=file_name2,
                    vcf_path=gl.get_path(vcf_file),
                    region=(chrom,region_start,region_end),
                    saveargs={"dpi":80,"facecolor":"white"}
                )
            else:
                mysumstats.plot_mqq(
                    save=file_name2,
                    region=(chrom,region_start,region_end),
                    saveargs={"dpi":80,"facecolor":"white"}
                )

    i +=1
        
#creating the parser
parser = argparse.ArgumentParser(description='Process multiple input paths')

#Adding the input path argument
parser.add_argument('--build',help='')
parser.add_argument('--chr_filter',help='')
parser.add_argument('--output_path', type=str, help='output path')
parser.add_argument('--input_paths', type=str, nargs='+', help='input paths')
parser.add_argument('--cut',type=int,help='cut value for plot_mqq')
parser.add_argument('--skip',type=int,default=0, help='skip value for plot_mqq')
parser.add_argument('--highlight',default=[],nargs="*", help='highlight value for plot_mqq')
parser.add_argument('--pinpoint',type=list, default=[],help='pinpoint value for plot_mqq')
parser.add_argument('--sig_level',default=5e-8,help='sig_level value for plot_mqq')
parser.add_argument('--pinpoint_color', type=str, default="red", help='pinpoint_color value for plot_mqq')
parser.add_argument('--highlight_color', type=str, default='#CB132D', help='highlight_color value for plot_mqq')
parser.add_argument('--anno',default="GENENAME",help='The variants to annotate will be selected automatically using a sliding window with windowsize=500kb')
parser.add_argument('--vcf_file',default=None,help="")

#parsing the arguments
args = parser.parse_args()

#processing the data
process_data(args.input_paths,args.output_path,args.cut,args.sig_level,args.skip,args.highlight,args.pinpoint, args.pinpoint_color,args.anno,args.highlight_color,args.build,args.chr_filter,args.vcf_file)
