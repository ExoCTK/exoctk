sudo yum -y install bzip2
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 700 ./Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
export PATH=/home/ec2-user/miniconda3/bin:$PATH
export EXOCTK_DATA=''
conda create --yes -n exoctk-aws python=3.6 git numpy flask pytest
conda init bash
source ~/.bashrc
conda activate exoctk-aws
git clone https://github.com/ExoCTK/exoctk.git
git clone https://github.com/spacetelescope/jwst_gtvt.git
cd jwst_gtvt
python setup.py develop
cd ../exoctk/
python setup.py develop