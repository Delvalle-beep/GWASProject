
                                  |G|W|A|S| |P|R|O|J|E|C|T|

<hr>
GWAS Project is a script that allows you to generate visualizations from summary statistics of a genome-wide study in a fully customized and automated way.
The purpose of its construction is to automate the researcher's time to carry out analyzes of a broad genome study, since it is a very complete and easy-to-use tool.
Its differential lies in having a direct connection with <a href="https://www.ensembl.org/index.html">Ensembl</a>, which facilitates the identification of significant variants of the study.

<h2>Quick Start</h2>

The mandatory arguments to run the program are "input_path","output_path" and "build"

Example:
```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' --build 'version of the genome being used in the genetic study'
```
if you want to run more than one file at the same time you just need to separete it with comma 

Example:
```
python GWASproject.py --input_path 'your/input/path/', 'your/input/path2/' ,'your/input/path3/' --output_path 'your/output/path' --build '19'
```
The console output should start like this:

![image.png](../tutorial-images/image.png)

<ul>
<li>Input files must be in "csv" or "tsv" format</li>
<li>The script will return a Manhattan plot for the GWAS summary statistics and a Regional plot for each significant variant in the study.</li>
</ul>

An example of the plots that can be generated with the script are:

![image1.png](../tutorial-images/image1.png)
![Regional1intro](../tutorial-images/Regional1intro.png)
![Regional2intro](../tutorial-images/Regional2intro.png)
![image3.png](../tutorial-images/image3.png)

The script is an optimization of the <a href="https://cloufield.github.io/gwaslab/">"GWASLAB"</a> library.
To get more personalized results some variables can be set, learn more in the next session.

<h2>Requirements</h2>

The packages required to run GWAS Project are:
<ul>
<li>matplotlib==3.6.3 or more updated version</li>
<li>gwaslab==3.3.18 or more updated version</li>
<li>pandas==1.4.4 or more updated version</li>
</ul>

This guide will help you get more unique and personalized views for your study.

<hr>

<h2>What is GWAS?</h2>
Genome-wide association studies (GWAS) are extremely important for identifying genetic polymorphisms that are associated with a specific outcome. Essentially, these tools scan specific chromosomes to identify which polymorphisms are associated with a given outcome. Consequently, there is an entire technology involved in identifying eligible polymorphisms that may be associated with a particular clinical study, for example.

The most common graph used for this type of study is called a Manhattan plot, and the degree of statistical significance required to identify a significant association between a polymorphism and a clinical outcome in GWAS is extremely conservative. That is, the P-value has to be very low, to the point of being in scientific notation, which is a requirement because it refers to another important concept in statistics. In the case of GWAS, since hundreds or thousands of polymorphisms are evaluated simultaneously, we need to correct the P-value using a technique called the "Bonferroni correction," where the P-value is multiplied according to each comparison. Therefore, the P-values need to be very low, in scientific notation, to be considered significant.

It is important to note that GWAS studies require a large sample size of thousands of individuals to be evaluated because we need high statistical power.

<h2>What is SNPs?</h2>
Single nucleotide polymorphisms (SNPs) are genetic variants that occur at a determined frequency in a population. This is the main difference between SNPs and simple mutations. SNPs occur at expected and determined frequencies, while mutations occur at unexpected frequencies. Furthermore, it is important to highlight that SNPs occur with a certain regularity in our genome. It is estimated that there is one SNP per thousand nucleotides, considering that our genome has 3 billion nucleotides, the estimate is that there are 3 million SNPs (although this can vary according to the study). More than 100 million SNPs that have some clinical or biological significance have been reported in the literature. It is important to note that many of these SNPs do not affect our organism or biology, as they can occur in regions of our DNA that will not affect gene expression. However, others may affect genes, generate important problems, metabolic disorders, and may be risk factors for diseases.

<h2>Study materials</h2>
![gwastutorialPrint](../tutorial-images/gwastutorialPrint.png)
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

![GWASLABorg](../tutorial-images/GWASLABorg.png)

Despite being originally written in Traditional Chinese, the site is well written and easily translated by public translation tools. I consider this the most complete and most explanatory material for beginners in the area.

<h1>Instructions for uset<h1>
<h2>Loading files</h2>
To load a file for reading in the script is simple, just inform the path to the file that you want to be read using the <b>input_path</b> command. 
For example:

```
python GWASproject.py --input_path 'your/input/path/' 
```
GWASProject was developed to be able to read multiple files with just one command, just inform them in this reading format:

```
python GWASproject.py --input_path 'your/input/path/', 'your/input/path2/' ,'your/input/path3/'
```
<b><i>Mandatorily, for the script to carry out the expected work, the input path must be defined, so that the study in question can be analyzed, an output path so that the generated files can be saved and shown, and finally, the version must be informed. of the genome on which the study was constructed.</i></b>

So that the script can read the columns of the loaded dataset, the column labels of this file must be informed. By default, GWASProject has a pre-definition of the most used labels, which are:

```         
['beta'],['p_value'],['snpid'],['variant_id'],['effect_allele'],['chromosome'],['other_allele'],
['standard_error'],['base_pair_location']          
```
![image4.png](../tutorial-images/image4.png)

These columns in turn are formatted within the script to the <a href="https://www.cog-genomics.org/plink/1.9/formats">PLINK</a> format.
Using the <a href="https://www.cog-genomics.org/plink/1.9/formats">PLINK</a> format, these columns change their label to:

```         
CHR=['chromosome']
POS=['base_pair_location']
rsID=['variant_id']
P=['p_value']
EA=['effect_allele']
NEA=['other_allele']
BETA=['beta']
SE=['standard_error']
SNPID=['snpid']
```     
![image5.png](../tutorial-images/image5.png)

<h2>File output</h2>
To define a path where the files generated after processing will be saved, you must define the  <b>output_path</b>. 
An application example is:
```
python GWASproject.py --input_path 'your/input/path/', 'your/input/path2/' ,'your/input/path3/' --output_path 'your/output/path' 
```
Unlike the previous topic, the <b>output_path</b> argument only receives a path where all the study files that were processed by the script will be saved.

It is important to point out that the files are saved according to their indexing in the order in which it was defined in the <b>input_path</b> argument.
For example, if you run the command:
```
python GWASproject.py --input_path 'first', 'second/' ,'third' 
```
In the folder defined in the <b>output_path</b>, it will be saved as:
```
“first1” , “second2”, “third3”
```

<h2>Genome build version</h2>
It refers to the specific representation of the human genome used as a reference to identify genetic variations associated with traits or diseases. 
The human genome is composed of a set of DNA sequences, and each sequence is identified by a specific position in the reference genome.
The genome build version is usually identified by a number, such as GRCh37 (also known as hg19) or GRCh38 (also known as hg38), which are the most commonly used versions. 
These designations correspond to different assemblies of the human genome, which are updated over time as new information is discovered about the structure of the genome.
Currently, GWASProject only supports the most common versions, which are hg19 and hg38, to set the version used in the summary statistics file, just use the <b>build</b> argument.

```
python GWASproject.py –build “19”
```

This argument, together with <b>input_path</b> and <b>output_path</b> form the mandatory arguments for the script to run correctly, without informing them the script will return an error.

```
python GWASproject.py --input_path 'your/input/path/', 'your/input/path2/' ,'your/input/path3/' --output_path 'your/output/path' –build “19”
```

<h2>Skip variants</h2>
The <b>skip</b> argument is used to skip a specified number of variants when generating a graph plotting <b>p-values</b> or <b>-log10(p)</b> values.
By supplying an integer value for the <b>skip</b> argument, the graph will be generated excluding the initial variants, skipping the specified number of variants. 
This can be useful when your dataset is very large and you want to reduce the density of variants displayed in the chart to improve readability.
For example, if you set <b>skip</b> to 10, then out of 10 variants, only one will be included in the graph. This can help reduce overlapping points on the graph, making it clearer and easier to interpret.
Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19” –skip 2
```
The Manhattan plot view for a <b>skip</b> value of 2 would be:

![image6.png](../tutorial-images/image6.png)

<h2>Cut in plotted values</h2>
The <b>cut</b> argument is used to define a cut value in plotting graphs of <b>p-values</b> or <b>-log10(p)</b> values. 
By supplying a numeric value for the <b>cut</b> argument, the plot will be limited to variants whose <b>p-values</b> or <b>-log10(p)</b> values are below this cut-off value. 
This can be useful for highlighting the most significant or relevant variants on the graph, filtering out the rest.
For example, if you set <b>cut</b> to 5, only variants with <b>p-values</b> below 5 (or <b>-log10(p)</b> values above 5) will be displayed on the graph. 
This allows you to focus on the strongest or most significant associations.

Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19” –cut 5
```

The Manhattan plot view for a cut value of 5 would be:

![image7.png](../tutorial-images/image7.png)

<h2>Significance level</h2>
The <b>sig_level</b> argument is used to define the level of statistical significance for plotting the <b>p-values</b> or <b>-log10(p)</b> values of the plotted graphs.
By providing a numeric value for the <b>sig_level</b> argument, the argument will be highlighted with a horizontal line representing the threshold of statistical significance. This helps identify which variants have <b>p-values</b> below this threshold and are considered statistically significant.

Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19” –sig_level 5e-6
```
For example, if you set <b>sig_level</b> to "5e-8", the graph will be plotted with a horizontal line representing a <b>p-value</b> of "5e-8". All variants with <b>p-values</b> below this threshold will be displayed above the line, indicating that they are statistically significant.

The Manhattan plot view for a sig_level value of 5e-8 would be:

![image1.png](../tutorial-images/image1.png)

<i><b>It is important to say that this argument only accepts float values, and that the default value for this argument is 5e-8.</b></i>

<h2>Highlight Argument</h2>
The <b>highlight</b> argument is used to highlight specific variants in the graphs plotted by the script. It allows you to provide a list of variants you want to highlight.

The <b>highlight_color</b> argument is used to set the highlight color for the variants specified in the <b>highlight</b> argument. You can provide a color in string format such as “red”, "#FF0000" (hex code) or "(1.0, 0.0, 0.0)" (RGB values).

When the variants specified in <b>highlight</b> are plotted on the chart, they are visually highlighted with the color specified in <b>highlight_color</b>, making it easier to identify these variants.

Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19” --highlight "10:69083:T:C" "10:94263:A:C" --highlight_color "#FF0000"
```

The Manhattan plot view for a value of <b>highligh</b> and <b>highlight_color</b> would be:

![image8.png](../tutorial-images/image8.png)

Make sure that the value chosen for the <b>highlight</b> argument is valid and is from the dataset being read.

<h2>Annotation</h2>
The <b>anno</b> argument is used to specify the annotation of variants on graphs plotted by <b>GWASProject</b>.

The <b>anno</b> refers to the column of the data that contains the annotation information for each variant.

This argument takes boolean, string or “GENENAME” values:

Boolean:
If the <b>anno</b> argument is set to <b>True</b>, the annotated variants will be automatically selected using a 500kb window and the chromosome:position format.
Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19”  --anno “True”
```

The Manhattan plot view for a “True” year value would be:

String:
the name of a column used for annotation will be used
Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19”  --anno “rsID”
```

The Manhattan plot view for an anno “rsID” value would be:

![image10.png](..tutorial-images/image10.png)

GENENAME:
This value is set to the default. Gene names are automatically annotated using <a href="https://github.com/openvax/pyensembl">pyensembl</a>.

Here is an example of how to use the argument:

```
python GWASproject.py --input_path 'your/input/path/' --output_path 'your/output/path' –build “19”  --anno “GENENAME”
```

The Manhattan plot view for an <b>anno</b> value “GENENAME” would be:
