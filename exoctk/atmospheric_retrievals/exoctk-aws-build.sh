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
echo '~~~~~ CREATING exoctk-aws CONDA ENVIRONMENT ~~~~~'
echo ''
conda create --yes -n exoctk-aws python=3.6 git numpy flask pytest
conda init bash
source ~/.bashrc
conda activate exoctk-aws

echo ''
echo '~~~~~ INSTALLING jwst_gtvt ~~~~~'
echo ''
git clone https://github.com/spacetelescope/jwst_gtvt.git
cd jwst_gtvt
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
python setup.py develop

echo ''
echo '~~~~~ INSTALLING ADDITIONAL LIBRARIES ~~~~~'
echo ''
pip install bibtexparser==1.1.0
pip install bokeh==1.3.1
pip install boto3==1.9.199
pip install corner==2.0.1
pip install h5py==2.8.0
pip install lmfit==0.9.13
pip install matplotlib==3.1.0
pip install pandas==0.25.0
pip install paramiko==2.4.2
pip install platon==3.1
pip install scp==0.13.2

echo ''
echo '~~~~~ THE ENVIRONMENT BEING USED ~~~~~'
echo ''
conda env export