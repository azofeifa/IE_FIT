#IE_FIT_1.1
##python dependencies
numpy, scipy, matplotlib
##module: convertAnnotation
The classify module takes in as input a specific annotation file format. "convertAnnotation" makes this for you however this is straightforward to make yourself.

$python NU_FIT/ convertAnnotation -i \<FILE.tsv> -j \<int int int int int> -o \<outFile.tsv>

-i: input file from USCS genome table, refseq annotation etc. must tab seperated
-j: provides the column number in 0-based coordinates for strand, chromosome, start coordinate, stop coordinate, geneName in that order
-o: output file name. The file outputted is strand[tab]chromosome[tab]start[tab]stop[tab]name/ID[newline]

##module: classify
Here we fit a mixture of a normal distribution followed by a uniform distribution to identify intiation and elongation events in genes. To estimate these parameters for each gene in your annotation file please look at the below command

$python NU_FIT/ classify -i \<genomeCoverageFile.bedgraph> -j \<annotationFormattedFile.tsv> -s \<strand, +/->  -o \<outputFileName.bed> \<...other non essential options...>

-i: \<path/to> input file, a genome coverage file must be formatted as chrom[tab]start[tab]stop[tab]coverage[newline]

-j: \<path/to> annotation file, the output from convertAnnotation or a file formatted as strand[tab]chromosome[tab]start[tab]stop[tab]name/ID[newline]

-o: \<path/to> ouputFile name \<any>

-s: \<string> strand either specify + or -; this tells "classify" which strand to pull from the annotation file and what orientation of the bedgraph file

-chr: \<string> (optional) you can specify a specific chromosome to run on, this string (i.e. chr1, chr2 ... ) must actually be in the genome coverage file and annotation file, default runs against all

-single: \<None> (optional) this will only take annotations in the annotation file that do not overlap another gene, i.e. single isoform genes, default behavior is off.

-np: \<int> (optional) number of processors, default is 8

-BIC:\<int int>, (optional, but recommended) first number specifies the max number of mixture models to run on each gene and the second int specifies the BIC criteria penality for increasing model selection. default is disabled and will run only one mixture model. 

-rt: \<int>, (optional) number of random iterations per each model fit, default is 1. 

-bin: \<int>, (optional) read coverage across annotation will be binned according to this resolution, default is 1000, higher the better

-merge: \<int>, (optional) this will take overlapping isoforms and take the min start and max stop of the isoforms as one single annotation

-int: \<int>, (optional) this takes values 0-20, if 0 then the first 5% of annotations will be calculated...this is useful for node parrellization i.e. job array in torque/qsub scripting

##module: concatenate
$python NU_FIT/ concatenate -i \<path/to/directory/of/IE_OUT files> -j \<path/to/out_file_name>

-i: \<path/to/> whatever the directory of IE_OUT files that you want to pull together and make into a single bed file

-j: \<path/to/> whatever you want to call the out file name, remember to append ".bed" to outfile name to load into IGV or UCSC genome browser

