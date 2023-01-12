# Toy model to test fluorescence subroutines within CORSIKA framework

Compile: `gfortran -o fluor src/test_fluorescence_corsika.F`

Execute: `./fluor` 

It should print out the number of fluorescence photons yielded at a given height and deposited energy defined in the source code.

The code is tested (compiled and executed) in a CI workflow (see [Actions](https://github.com/morcuended/test-fluorescence/actions/workflows/ci.yml) tab) everytime a change is introduced.

[![CI](https://github.com/morcuended/test-fluorescence/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/morcuended/test-fluorescence/actions/workflows/ci.yml)
