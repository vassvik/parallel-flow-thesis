parallel-flow-thesis
====================

Source code used in master thesis simulations

Compile using

g++ *.cpp *.hpp -std=c++11 -O3

run with

cat config.txt | ./a.out >> log &

or using overriding command line arguments

cat config.txt | ./a.out -seed 1 -pressDiff 3.0 >> log &

Output is redirected to log-file.
