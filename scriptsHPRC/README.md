#create conda env (once)
conda init
conda env create -f $HP_SDIR/scripts/install_prerequisites.yaml  

# activate
conda activate mitohpc2

# init (set the HP_ variables)
$HP_SDIR/init.sh
$HP_SDIR/init.lr.sh
$HP_SDIR/init.ont.sh

# check all vars are set
$HP_SDIR/checkInstall.sh                       
$HP_SDIR/checkInstall.ref.sh                

############################################


mkdir Illumina
cd Illumina

# download one Illumina sample
head -1 $HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh 
  download_aws_illumina.sh HG00097	s3://human-pangenomics/submissions/59C50DDF-5FAF-4841-AC3E-6C02D636C57F--Y4_1000G_DATA/HG00097.final.cram
head -1 $HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh | bash

# download all Illumina samples
$HP_SDIR/../scriptsHPRC/download_aws_illumina.all.sh

