#!/bin/bash
# Script to setup environment
sudo apt update
sudo apt -y install make git clang++-10 libgmp-dev g++ parallel
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-10 1000
git submodule init
git submodule update
make -C ./include/ate-pairing

