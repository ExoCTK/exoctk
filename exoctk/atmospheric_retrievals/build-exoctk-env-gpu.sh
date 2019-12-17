echo ''
echo '~~~~~ INSTALLING NVIDIA GPU DRIVER ~~~~~'
echo ''
sudo yum -y update
sudo yum -y install wget nano elfutils-libelf-devel
sudo yum -y groupinstall "Development Tools"
sudo sed -i 's/crashkernel=auto"/crashkernel=auto nouveau.modeset=0"/g' /etc/default/grub
sudo grub2-mkconfig -o /boot/grub2/grub.cfg
sudo touch /etc/modprobe.d/blacklist.conf
sudo chmod 777 /etc/modprobe.d/blacklist.conf
sudo echo 'blacklist nouveau' > /etc/modprobe.d/blacklist.conf
sudo mv /boot/initramfs-$(uname -r).img /boot/initramfs-$(uname -r)-nouveau.img
sudo dracut /boot/initramfs-$(uname -r).img $(uname -r)
sudo reboot
sudo systemctl isolate multi-user.target
wget http://us.download.nvidia.com/XFree86/Linux-x86_64/430.40/NVIDIA-Linux-x86_64-430.40.run
sudo sh NVIDIA-Linux-x86_64-430.40.run
sudo reboot

echo ''
echo '~~~~~ INSTALLING CUDA TOOLKIT ~~~~~'
echo ''
wget http://developer.download.nvidia.com/compute/cuda/10.1/Prod/local_installers/cuda_10.1.243_418.87.00_linux.run
sudo sh cuda_10.1.243_418.87.00_linux.run
export PATH=$PATH:/usr/local/cuda-10.1/bin
export LD_LIBRARY_PATH=/usr/local/cuda-10.1/lib64

echo ''
echo '~~~~~ INSTALLING ANACONDA ~~~~~'
echo ''
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 700 ./Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

echo ''
echo '~~~~~ SETTING ENVIRONMENT VARIABLES ~~~~~'
echo ''
export PATH=/home/ec2-user/miniconda3/bin:$PATH
export EXOCTK_DATA=''

echo ''
echo '~~~~~ CREATING base CONDA ENVIRONMENT ~~~~~'
echo ''
conda create --yes -n exoctk-3.6 python=3.6 git numpy flask pytest
conda init bash
source ~/.bashrc
conda activate exoctk-3.6

echo ''
echo '~~~~~ INSTALLING exoctk ~~~~~'
echo ''
git clone https://github.com/ExoCTK/exoctk.git
cd exoctk/
conda env update -f env/environment-3.6.yml
conda init bash
source ~/.bashrc
conda activate exoctk-3.6
python setup.py develop
cd ../

echo ''
echo '~~~~~ INSTALLING jwst_gtvt ~~~~~'
echo ''
rm -fr /home/ec2-user/miniconda3/envs/exoctk-3.6/lib/python3.6/site-packages/jwst_gtvt
git clone https://github.com/spacetelescope/jwst_gtvt.git
cd jwst_gtvt
git checkout cd6bc76f66f478eafbcc71834d3e735c73e03ed5
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
echo '~~~~~ INSTALLING cudamat ~~~~~'
echo ''
git clone https://github.com/cudamat/cudamat.git
sudo sed -i "s/-O',/-O3',/g" /home/ec2-user/cudamat/setup.py
cd cudamat
python setup.py develop
cd ../

echo ''
echo '~~~~~ INSTALLING gnumpy ~~~~~'
echo ''
git clone https://github.com/ExoCTK/gnumpy3.git
cd gnumpy3
python setup.py develop
cd ../

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