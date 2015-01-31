#PBS -N IE_FIT
#PBS -l nodes=1:ppn=64
#PBS -t 1-20
module load numpy_1.8.1
module load scipy_0.14.0 
module load python_2.7.3

root=/Users/azofeifa/
outFilePath=${root}outFiles_IE_FIT/
inFilePath=${root}inFiles_IE_FIT/
pathTosrc=${root}IE_FIT/NU_FIT/

BG_FILE=DMSO2_3.sorted.fiveprime.pos.BedGraph
annot_FILE=RefSeq.NU.tsv
strand=+
chrom=chr1
single=1
test=0
rt=5
BIC_MAX=5
BIC_PEN=$PBS_ARRAYID
outFile=${chrom}_${BIC_MAX}_${BIC_PEN}_${rt}_IE_OUT.bed
np=64
if [ $test == "1" ]
then
    python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom -v -BIC $BIC_MAX $BIC_PEN -np $np -t
else
    if [ $single == "1" ]
    then
	python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom  -BIC $BIC_MAX $BIC_PEN -single -np $np -rt $rt
    else
	python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom  -BIC $BIC_MAX $BIC_PEN -np $np -rt $rt -merge
    fi
fi



