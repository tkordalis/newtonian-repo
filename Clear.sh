#!/bin/bash

rm -rf *.plt
rm -rf pp0*
rm -rf *.mod
rm -rf *.bck
rm -rf Boundaries
rm -rf solution
#rm -rf tecplot
rm -rf tags
rm tags
rm -rf nohup.out
rm -rf .read_unv.py
rm -rf *.PLT
rm -rf *.DAT
rm -rf *.dat
rm -f  ./BubbleVariables/*
rm -f  ./BubbleVariables/Bubble1/*
rm -f  ./BubbleVariables/Bubble2/*
rm -f  *_genmod.f90
rm -rf  .UnBounded.unv/
rm -rf  .time_*.unv/
ctags -R .
