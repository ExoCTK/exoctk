echo ''
echo '~~~~~ INSTALLING ANACONDA ~~~~~'
echo ''
sudo yum -y install bzip2
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 700 ./Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

echo ''
echo '~~~~~ SETTING ENVIRONMENT VARIABLES ~~~~~'
echo ''
export PATH=/home/ec2-user/miniconda3/bin:$PATH
echo 'Set $PATH'
export EXOCTK_DATA=''
echo 'Set $EXOCTK_DATA'

echo ''
echo '~~~~~ CREATING base CONDA ENVIRONMENT ~~~~~'
echo ''
conda create --yes -n exoctk-3.6 python=3.6 git numpy flask pytest
conda init bash
source ~/.bashrc
conda activate exoctk-3.6

echo ''
echo '~~~~~ INSTALLING jwst_gtvt ~~~~~'
echo ''
git clone https://github.com/spacetelescope/jwst_gtvt.git
cd jwst_gtvt
git checkout cd6bc76f66f478eafbcc71834d3e735c73e03ed5
python setup.py develop
cd ../

echo ''
echo '~~~~~ INSTALLING exoctk ~~~~~'
echo ''
git clone https://github.com/ExoCTK/exoctk.git
cd exoctk/
git remote add bourque https://github.com/bourque/exoctk.git
git fetch bourque
git checkout -b implement-aws bourque/implement-aws
conda env update -f env/environment-3.6.yml
conda init bash
source ~/.bashrc
conda activate exoctk-3.6
python setup.py develop
cd ../

echo ''
echo '~~~~~ INSTALLING ADDITIONAL LIBRARIES ~~~~~'
echo ''
pip install bibtexparser==1.1.0
pip install corner==2.0.1
pip install lmfit==0.9.13
pip install platon==3.1

echo ''
echo '~~~~~ BUILDING EXOCTK_DATA DIRECTORY ~~~~~'
echo ''
mkdir exoctk_data/
mkdir exoctk_data/exoctk_contam/
mkdir exoctk_data/exoctk_log/
mkdir exoctk_data/fortney/
mkdir exoctk_data/generic/
mkdir exoctk_data/groups_integrations/
mkdir exoctk_data/modelgrid/

echo ''
echo '~~~~~ THE ENVIRONMENT BEING USED ~~~~~'
echo ''
conda env export