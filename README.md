# GWASproject


## Getting started

To run the program you must set two arguments "input_path" and "output_path"

Example:
```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path'
```
if you want to run more than one file at the same time you just need to separete it with comma 

Example:
```
python GWASproject.py --input_path 'your/input/path/', 'your/input/path2/' ,'your/input/path3/' --output_path 'your/output/path'
```

The script is an optimization of the <a href="https://cloufield.github.io/gwaslab/">GWASLAB"</a> library.
To obtain more personalized results in the generated visualizations, some variables can be set.

#Manhattan and QQ Plot (plot_mqq)

For the QQ and Manhattan plots the default values ​​are the following:
```
mysumstats.plot_mqq(
          anno=True,
          cut=0,
          skip=0,
          sig_level=5e-8,
          highlight = [],
          highlight_color="#CB132D",
          pinpoint=[],
          pinpoint_color ="red",
          marker_size=(5,25),
          build="19",
          region_mqq=None,
          saveargs={"dpi":400,"facecolor":"white"}
          )
```
<ul>
<li><b>Skip:</b> sometimes it is not necessary to plot all variants, we can skip the insignicant variants . For example, we can exclude varints with -log10p lower than 3 from the plot by specifying skip=3</li>
Use this command to set:
</br>
```
--skip "0"
```
</br>
<li><b>Cut:</b>loci with extremly large -log10(P) value are very likely to dwarf other significant loci , so we want to scale down the extrame loci from a certain threshold.</li>
Use this command to set:
</br>
```
--cut "0"
```
</br>
<li><b>Sig_level:</b>genome-wide significance threshold. Specify the P value threshold.</li>
</br>
Use this command to set:
```
--sig_level "5e-8"
```
</br>
<li><b>Anno:</b>When it is set "True", the variants to annotate will be selected automatically using a sliding window with windowsize=500kb.</li>
</br>
Use this command to set:
```
--anno True
```
</br>
<li><b>Highlight:</b>A "type=list" element to specify the variants of loci for highlighting.</li>
</br>
Use this command to set:
```
--highlight ["rs12509595","19:15040733:T:C"]
```
</br>
<li><b>Highlight_color:</b>specify the color ussed for highlighting, RGB format.</li>
</br>
Use this command to set:
```
--highlight_color "#CB132D"
```
</br>
<li><b>Pinpoint:</b>a list of SNPIDs.</li>
</br>
Use this command to set:
```
--pinpoint ["rs7989823"]
```
</br>
<li><b>Pinpoint_color:</b>color for pinpoint.</li>
</br>
Use this command to set:
```
--pinpoint_color "red"
```
</br>
<li><b>Build:</b>genome build version "19" or "38".</li>
</br>
<li><b>Region_mqq:</b>A "type=tuple" element to define wich region to plot. Example: region=(7,156538803,157538803)</li>
</br>
Use this command to set:
```
--region_mqq (7,,156538803,157538803) 
```
</br>
<li><b>marker_sizer:</b>It define the marker size,is defined by two integers, an example of use is:marker_size=(5,10)</li>
</br>
Use this command to set:
```
--marker_size (5,25)
```
</br>
</ul>

#Regional Plot

For the Regional Plot  the default values ​​are the following:
```
mysumstats.plot_mqq(
          anno=True,
          cut=0,
          skip=0,
          sig_level=5e-8,
          highlight = [],
          highlight_color="#CB132D",
          pinpoint=[],
          pinpoint_color ="red",
          region_regPlot=None,
          vcf_path=None,
          vcf_chr_dict=None,
          mode="r"
          )
```
Many of the arguments are also used in the Manhattan and QQ plots and have the same definition, the definition of additional arguments will be added below.

<ul>
<li><b>Region_regPlot:</b> sometimes it is not necessary to plot all variants, we can skip the insignicant variants . For example, we can exclude varints with -log10p lower than 3 from the plot by specifying skip=3. This argument is used to uniquely define the region plotted on regional plots only.</li>
Use this command to set:
</br>
```
--region_regPlot (7,,156538803,157538803)
```
</br>
<li><b>Vcf_path:</b>If you want to load a reference genotype just add the path</li>
Use this command to set:
</br>
```
--vcf_path 'C:/This/Is/My/Path'
```
</br>
<li><b>Vcf_chr_dict:</b>This argument is defined to select one or a portion of chromosomes to be used in the reference genotype. This already has a built-in get() function that filters the chromosome you would like to reference.</li>
</br>
Use this command to set:
```
--vcf_chr_dict '6'
```
</br>
</ul>
