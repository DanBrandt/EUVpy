#!/bin/sh

gfortran -c -ffixed-line-length-132 HEUVAC-Driver.for
gfortran -c -ffixed-line-length-132 HEUVAC.for
gfortran -o HEUVAC.exe HEUVAC-Driver.o HEUVAC.o
