#!/bin/bash
cd ./src &&
f2py -m -c mstarsim mstarsim.f95 --opt=-O3 &&
mv mstarsim.so ../bin
