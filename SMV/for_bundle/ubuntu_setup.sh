#!/bin/bash
USERNAME=`whoami`

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

INSTALL openssh-server
INSTALL samba

# smbpasswd -a gforney
# add following to end of /etc/samba/smb.conf
#[gforney]
#    path = /home/gforney
#    browsable = yes
#    guest ok = yes
#    read only = no
#    create mask = 0755

# sudo restart smbd
# sudo restart nmbd

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
INSTALL g++multilib >> setup_ubuntu.txt
INSTALL freeglut3
INSTALL freeglut3-dev
INSTALL texlive-latex-base
INSTALL texlive-latex-extra
INSTALL texlive-science
INSTALL texlive-fonts-recommended
INSTALL xvfb
INSTALL msmtp
INSTALL heirloom-mailx

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
EOF

cat << EOF > ~/.mailrc_ORIG
# --- ~/.mailrc ---
set sendmail="/usr/bin/msmtp"
set message-sendmail-extra-arguments="-a gmail"
EOF
 
#account default : gmail

echo adding grive repository
echo adding grive repository >> setup_ubuntu.txt
echo >> setup_ubuntu.txt
add-apt-repository -y ppa:thefanclub/grive-tools >> setup_ubuntu.txt
apt-get -y update >> setup_ubuntu.txt
INSTALL grive-tools

echo setup complete
