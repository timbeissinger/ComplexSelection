#BSUB -J Run_50
#BSUB -o bsub_log/Run_50.o%J
#BSUB -e bsub_log/Run_50.e%J
#BSUB -q normal
#BSUB -n 16
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -u beissingert@missouri.edu

### Set variables
qmsimDir=/home/beissingert/QMSim_Linux/
paramFile=qtl_50.prm
dir=r_qtl_50

### Run QMSim
$qmsimDir/QMSim $paramFile -o

### Trim genotypes file to 1000 individuals
for i in $(ls $dir/DriftPop_mrk_*.txt)
do
    head -n 1001 $i > $i-head
    rm $i
done

for i in $(ls $dir/SelectPop_mrk_*.txt)
do
    head -n 1001 $i > $i-head
    rm $i
done



