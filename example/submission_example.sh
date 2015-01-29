

########## first the high level parameters about where the scripts and the bundle are located
RNASEQPIPBASE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline  ##this probably needs to be updated
RNASEQBUNDLE=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle  ##this probably should stay



export RNASEQPIPBASE=$RNASEQPIPBASE
export RNASEQBUNDLE=$RNASEQBUNDLE
pipeline=${RNASEQPIPBASE}/RNAseq_pipeline_v7.sh

echo "Main pipeline script is hereL $pipeline"

###################
##All the parameters to run the pipeline are below

submit=yes
force=no

tophat=no
sampleQC=no
dexseqcounts=no
runCufflinks=no
miso=no

prepareCounts=no
Rdeseq=no
Rdexseq=yes

RpathwayGO=no
RtopGO=no

code=Zanda_AD_Tc1J20  ###identifier for the run

iFolder=/SAN/biomed/biomed14/vyp-scratch/Zanda_AD_Tc1J20_RNASeq/fastq/  ##input folder that contains the fastq files
oFolder=/scratch2/vyp-scratch2/IoN_RNASeq/Frances/processed  ### output folder that will contain the output data

species=tc1_mouse ## choice of species used

mainscript="cluster/submission/de_$code.sh"  ##main script, we will apply qsub to this
dataframe=data/RNASeq_AD_Tc1J20.tab  ##main and only support file, describing conditions and covariates

bash $pipeline --tophat ${tophat} --dexseqcounts ${dexseqcounts} --runCufflinks ${runCufflinks} --sampleQC ${sampleQC} --mainscript $mainscript --iFolder ${iFolder} --oFolder ${oFolder} --dataframe $dataframe --code $code --prepareCounts ${prepareCounts} --Rdexseq ${Rdexseq} --Rdeseq ${Rdeseq} --RpathwayGO ${RpathwayGO} --RtopGO ${RtopGO} --submit ${submit} --species $species --force $force 

echo $mainscript
