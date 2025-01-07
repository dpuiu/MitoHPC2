#!/usr/bin/env bash
set -x

##############################################################################################################

# Program that checks/installs packages needed by bwa, samtools, bedtools ...
# To be run using sudo; sudo privileges required !!!

##############################################################################################################

command -v apt-get
if [ "$?" == 0 ] ; then
  apt-get -y update # && apt-get upgrade
  apt-get install -y git wget openjdk-17-jdk zlib1g libz-dev libncurses5-dev libbz2-dev pkg-config liblzma-dev build-essential unzip  parallel # make gcc python
  apt-get install -y libcurl4-gnutls-dev libssl-dev zlib1g-dev
  apt-get install -y python-is-python3
  apt-get install -y gfortran libreadline-dev libpcre2-dev  # for R
fi

command -v dnf
if [ "$?" == 0 ] ; then
  dnf -y update
  dnf install -y which nano git wget java-17-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel unzip perl perl-Data-Dumper perl-ExtUtils-MakeMaker perl-Test-Simple python3 python3-pip make libcurl-devel openssl-devel
  alternatives --install /usr/bin/python python /usr/bin/python3 60
fi

command -v yum
if [ "$?" == 0 ] ; then
  yum -y update
  #yum install -y which nano git wget java-1.8.0-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel  unzip perl perl-Data-Dumper  perl-ExtUtils-MakeMaker perl-Test-Simple python parallel
  yum install -y which nano git wget java-17-openjdk bzip2 gcc gcc-c++ zlib-devel ncurses-devel bzip2-devel xz-devel  unzip perl perl-Data-Dumper  perl-ExtUtils-MakeMaker perl-Test-Simple python3 python3-pip make libcurl-devel openssl-devel # removed parallel python
  alternatives --install /usr/bin/python python /usr/bin/python3 60
fi
