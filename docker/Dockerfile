FROM ubuntu:22.04

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
RUN sudo apt-get install -yq make cmake gcc g++ gfortran libmpich-dev libopenblas-dev
RUN sudo apt-get install -yq vim
RUN sudo apt-get install -yq python3 python3-pip
# RUN sudo apt-get install -yq curl
# RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
# RUN sudo python2 get-pip.py
RUN sudo apt-get install -yq git-lfs
RUN sudo apt-get install -yq valgrind
RUN sudo apt-get install -yq hdf5-tools

RUN echo "numpy" >> requirements.txt
RUN echo "scipy" >> requirements.txt
# RUN cp requirements.txt requirements.2.txt
# RUN sudo pip install -r ./requirements.2.txt

RUN echo "argparse" >> requirements.txt
RUN echo "tables" >> requirements.txt
RUN echo "pandas" >> requirements.txt
RUN echo "PyYAML" >> requirements.txt
RUN echo "h5py" >> requirements.txt
RUN sudo pip3 install --upgrade pip
RUN sudo pip3 install -r ./requirements.txt

# Build mpi4py and petsc4py from source against the system MPICH.
# Default PyPI wheels target OpenMPI; ABI mismatch with libmpich-dev would
# silently corrupt MPI_Comm handles. petsc4py pulls in `petsc` as a sibling
# pip package and compiles PETSc itself from source.
RUN sudo pip3 install --no-binary=mpi4py mpi4py
RUN sudo pip3 install --no-binary=petsc,petsc4py petsc petsc4py

# Smoke-test the install: fail the build if PETSc + mpi4py can't initialize.
RUN python3 -c "from mpi4py import MPI; from petsc4py import PETSc; PETSc.Sys.Print('petsc4py/mpi4py ok: PETSc %d.%d.%d' % PETSc.Sys.getVersion())"

RUN sudo apt-get clean -q
