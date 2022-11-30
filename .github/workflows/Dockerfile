FROM ubuntu:18.04

ENV USERNAME=test

# install sudo
RUN apt-get -yq update && apt-get -yq install sudo

# create and switch to a user
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN useradd --no-log-init --home-dir /home/$USERNAME --create-home --shell /bin/bash $USERNAME
RUN adduser $USERNAME sudo
USER $USERNAME
WORKDIR /home/$USERNAME

# install packages
RUN sudo apt-get install -yq git
RUN sudo apt-get install --no-install-recommends -yq make cmake gcc gfortran openmpi-bin openmpi-common libopenmpi-dev
RUN sudo apt-get install -yq vim
RUN sudo apt-get install -yq python python3
RUN sudo apt-get install -yq python3-pip
RUN sudo apt-get install -yq curl
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN sudo python2 get-pip.py

RUN echo "numpy" >> requirements.txt
RUN echo "scipy" >> requirements.txt
RUN cp requirements.txt requirements.2.txt

RUN echo "argparse" >> requirements.txt
RUN echo "pandas" >> requirements.txt
RUN echo "PyYAML" >> requirements.txt

RUN sudo pip install -r ./requirements.2.txt
RUN sudo pip3 install -r ./requirements.txt
#RUN sudo pip install numpy
#RUN sudo pip install scipy
#RUN sudo pip3 install numpy
#RUN sudo pip3 install scipy
#RUN sudo pip3 install argparse
#RUN sudo pip3 install pandas
RUN sudo apt-get install -yq valgrind
RUN sudo apt-get install -yq git-lfs
RUN sudo apt-get clean -q

# get modern-fortran code
#RUN git clone https://github.com/modern-fortran/tsunami

