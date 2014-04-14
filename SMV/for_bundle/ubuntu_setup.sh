#!/bin/bash
# run this script using sudo

# set user name of person "owning" cfast and FDS-SMV repositories
USERNAME=gforney


INSTALL ()
{
  package=$1
  echo installing $package
  echo installing $package >> setup_ubuntu.txt
  echo >> setup_ubuntu.txt
  apt-get -y install $package >> setup_ubuntu.txt
  echo ------------------------------- >> setup_ubuntu.txt
}

echo > setup_ubuntu.txt
apt-get update >> setup_ubuntu.txt

INSTALL openssh-server
INSTALL samba

echo setting up samba
smbpasswd -a $USERNAME

cat << EOF >> /etc/samba.conf
[$USERNAME]
    path = /home/$USERNAME
    browsable = yes
    guest ok = yes
    read only = no
    create mask = 0755
EOF
restart smbd
restart nmbd

INSTALL subversion

echo installing cfast repository
echo installing cfast repository >> setup_ubuntu.txt
svn co https://cfast.googlecode.com/svn/trunk/cfast/trunk ~/cfast --username $USERNAME >> setup_ubuntu.txt
chown -R $USERNAME.$USERNAME ~/cfast
echo >> setup_ubuntu.txt
echo ------------------------------- >> setup_ubuntu.txt

echo installing FDS-SMV repository
echo installing FDS-SMV repository >> setup_ubuntu.txt
svn co https://fds-smv.googlecode.com/svn/trunk/FDS/trunk ~/FDS-SMV --username $USERNAME >> setup_ubuntu.txt
chown -R $USERNAME.$USERNAME ~/FDS-SMV
echo >> setup_ubuntu.txt
echo ------------------------------- >> setup_ubuntu.txt

INSTALL libxmu-dev
INSTALL libxi-dev
INSTALL g++-multilib
INSTALL freeglut3
INSTALL freeglut3-dev
INSTALL texlive-latex-base
INSTALL texlive-latex-extra
INSTALL texlive-science
INSTALL texlive-fonts-recommended
INSTALL xvfb
INSTALL msmtp
INSTALL heirloom-mailx
INSTALL gfortran
INSTALL mjpegtools

echo copy ~/.msmtprc_ORIG to ~/.msmtprc and edit
cat << EOF > ~/.msmtprc_ORIG
#---- ~/.msmtprc -----
#Gmail account
defaults
logfile ~/msmtp.log
 
account gmail
auth on
host smtp.gmail.com
from user@gmail.com
auth on
tls on
# tls_trust_file /etc/ssl/cert.pem
tls_trust_file /usr/share/ca-certificates/mozilla/Equifax_Secure_CA.crt
user user@gmail.com
password userpasswd
port 587

account default: gmail
EOF

echo copy ~/.mailrc_ORIG to ~/.mailrc
cat << EOF > ~/.mailrc_ORIG
# --- ~/.mailrc ---
set sendmail="/usr/bin/msmtp"
set message-sendmail-extra-arguments="-a gmail"
EOF

echo adding grive repository
echo adding grive repository >> setup_ubuntu.txt
echo >> setup_ubuntu.txt
add-apt-repository -y ppa:thefanclub/grive-tools >> setup_ubuntu.txt
apt-get update >> setup_ubuntu.txt
INSTALL grive
INSTALL grive-tools
INSTALL tcsh
INSTALL libnuma-dev
INSTALL gparted
INSTALL mesa-utils
INSTALL apache2
INSTALL gpicview

# adding NVIDIA video drivers (by hand)

# first time

#sudo apt-get install nvidia-current

# next and succeeding times

#sudo apt-get dist-upgrade

# install openmpi libraries

wget https://bintray.com/nist-fire-research/releases/Open_MPI/view/files/openmpi_1.6.5_linux_32.tar.gz
wget https://bintray.com/nist-fire-research/releases/Open_MPI/view/files/openmpi_1.6.5_linux_64.tar.gz
mkdir /shared
cd /shared
tar xvf ~/openmpi_1.6.5_linux_32.tar.gz
tar xvf ~/openmpi_1.6.5_linux_64.tar.gz

echo setup complete
