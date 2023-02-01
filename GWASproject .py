
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import gwaslab as gl
import pandas as pd
import sys
import os

user_input = input("Enter the path of your file: ")
user_output = input("Enter the output path of your file: ")  
input_file_name = input ("Enter the output file name: ")

assert os.path.exists(user_input), "I did not find the file at, "+str(user_input)
assert os.path.exists(user_output), "This output path did not exist "+str(user_output)
f = pd.read_csv(user_input, sep="\t",nrows=100000)
print("Hooray we found your file!")

mysumstats = gl.Sumstats(user_input,
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


user_output2 = user_output + '/' + input_file_name + '.png'
mysumstats.plot_mqq(save= user_output2,saveargs={"dpi":400,"facecolor":"white"})
mysumstats.get_lead()
user_output3 = user_output + '/' + input_file_name + '2' + '.png'
mysumstats.plot_mqq(save= user_output3,region=(7,156538803,157538803))
user_output4 = user_output + '/' + input_file_name + '3' + '.png'
mysumstats.plot_mqq(save= user_output4, mode="r", region=(7,156538803,157538803),region_grid=True, gtf_path="ensembl")

print("The script is over!! Go to the output folder and check the results!!")


