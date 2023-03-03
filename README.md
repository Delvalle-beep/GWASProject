
+-+-+-+-+ +-+-+-+-+-+-+-+
|G|W|A|S| |P|R|O|J|E|C|T|
+-+-+-+-+ +-+-+-+-+-+-+-+

GWAS Project is a script that allows you to generate visualizations from summary statistics of a genome-wide study in a fully customized and automated way.
The purpose of its construction is to automate the researcher's time to carry out analyzes of a broad genome study, since it is a very complete and easy-to-use tool.

<h2>Quick Start</h2>

The mandatory arguments to run the program are "input_path" and "output_path"

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
To get more personalized results some variables can be set. Learn more in the #GettingStarted session

<h2>Getting Started</h2>

<b>A practical tool for working with GWAS summary statistics</b>
GWAS Project is a script that allows you to generate visualizations from summary statistics of a genome-wide study in a fully customized and automated way.
The purpose of its construction is to automate the researcher's time to carry out analyzes of a broad genome study, since it is a very complete and easy-to-use tool.

The packages required to run GWAS Project are:
<ul>
<li>matplotlib==3.6.3 or more updated version</li>
<li>gwaslab==3.3.18 or more updated version</li>
<li>pandas==1.4.4 or more updated version</li>
</ul>

This guide will help you get more unique and personalized views for your study.

--------------------------------------------------------------------------------------------------

<h3>What is GWAS?</h3>
Genome-wide association studies (GWAS) are extremely important for identifying genetic polymorphisms that are associated with a specific outcome. Essentially, these tools scan specific chromosomes to identify which polymorphisms are associated with a given outcome. Consequently, there is an entire technology involved in identifying eligible polymorphisms that may be associated with a particular clinical study, for example.

The most common graph used for this type of study is called a Manhattan plot, and the degree of statistical significance required to identify a significant association between a polymorphism and a clinical outcome in GWAS is extremely conservative. That is, the P-value has to be very low, to the point of being in scientific notation, which is a requirement because it refers to another important concept in statistics. In the case of GWAS, since hundreds or thousands of polymorphisms are evaluated simultaneously, we need to correct the P-value using a technique called the "Bonferroni correction," where the P-value is multiplied according to each comparison. Therefore, the P-values need to be very low, in scientific notation, to be considered significant.

It is important to note that GWAS studies require a large sample size of thousands of individuals to be evaluated because we need high statistical power.

<h3>What is SNPs?</h3>
Single nucleotide polymorphisms (SNPs) are genetic variants that occur at a determined frequency in a population. This is the main difference between SNPs and simple mutations. SNPs occur at expected and determined frequencies, while mutations occur at unexpected frequencies. Furthermore, it is important to highlight that SNPs occur with a certain regularity in our genome. It is estimated that there is one SNP per thousand nucleotides, considering that our genome has 3 billion nucleotides, the estimate is that there are 3 million SNPs (although this can vary according to the study). More than 100 million SNPs that have some clinical or biological significance have been reported in the literature. It is important to note that many of these SNPs do not affect our organism or biology, as they can occur in regions of our DNA that will not affect gene expression. However, others may affect genes, generate important problems, metabolic disorders, and may be risk factors for diseases.

<h4>Study materials</h4>
![gwastutorialPrint](../uploads/c2d29e8d675ff08372e80a5780124686/gwastutorialPrint.png)
If population genetics or genome-wide association studies are new to you, I recommend visiting this <a href="https://cloufield.github.io/GWASTutorial/">website</a> for more in-depth content.

This tutorial is provided by the Kamaya Laboratory at the University of Tokyo. It is primarily intended for beginners in bioinformatics/statistical genetics. It covers the following parts:

<ul>
<li>Command line (mostly linux, a small amount of R/Python/JupyterNotebook/Github, etc.)</li>
<li>Data processing and quality control before GWAS</li>
<li>GWAS and results visualization</li>
<li>Downstream analysis after GWAS</li>
<li>GWAS Related Topics</li>
</ul>

To delve even further into statistical and computational concepts, I strongly recommend accessing this <a href="https://gwaslab.org/">site</a>.

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
