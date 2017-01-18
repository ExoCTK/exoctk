#! /bin/csh -fv
gfortran cea2.f -o cea2.x -w

