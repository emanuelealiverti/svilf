#!/bin/bash

#1] attempt to build and compile svilf in the local directory
# TODO: add a check, and compile only if svilf does not exists OR gives some error due to compilation on a wrong architecture
cwd=$(pwd)
cd ../
R -e 'if(!require("svilf", lib.loc = "./lib")) system("./R_localbuild svilf")'

#2] Compile stan competitor
cd $cwd/00_COMPILATION
R -e 'if(!file.exists("factor.RData")) source("./compile_stan.r")'

#3] check for the packages implementing latent distance
R -e 'if(!require("lvm4net")) install.packages("lvm4net")'
R -e 'if(!require("VBLPCM")) install.packages("VBLPCM")'

