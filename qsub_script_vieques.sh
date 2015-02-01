#PBS -N IE_FIT
#PBS -l nodes=1:ppn=32
#PBS -e /Users/azofeifa/qsub_errors/
#PBS -o /Users/azofeifa/qsub_stdo/
#PBS -t 1-500%20
module load numpy_1.6.1
module load scipy_0.12.0
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
rt=10
BIC_MAX=10
BIC_PEN=$PBS_ARRAYID
BIN=500
outFile=${chrom}_${BIC_MAX}_${BIC_PEN}_${rt}_${BIN}_IE_OUT.bed
np=32
if [ $test == "1" ]
then
    python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom -v -BIC $BIC_MAX $BIC_PEN -np $np -t
else
    if [ $single == "1" ]
    then
	python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom  -BIC $BIC_MAX $BIC_PEN -single -np $np -rt $rt -bin $BIN
    else
	python $pathTosrc classify -i ${inFilePath}$BG_FILE -j ${inFilePath}$annot_FILE -o ${outFilePath}$outFile -s $strand -chr $chrom  -BIC $BIC_MAX $BIC_PEN -np $np -rt $rt -merge
    fi
fi



