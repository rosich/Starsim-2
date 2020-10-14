#!/bin/bash
cd ./src &&
f2py2.7 -m -c mstarsim mstarsim.f95 --opt=-O3 &&
mv mstarsim*.so ../bin/mstarsim.so
