parallel-flow-thesis
====================

Source code used in master thesis simulations

Compile using

<code>
g++ *.cpp *.hpp -std=c++11 -O3
</code>

Modify config.txt accordingly, then run with

<code>
cat config.txt | ./a.out >> log &
</code>

or using overriding command line arguments

<code>
cat config.txt | ./a.out -seed 1 -pressDiff 3.0 >> log &
</code>

Output is redirected to the file 'log'.

You can continually read the output using tail, which will update when new output has been written:

<code>
tail -f log
</code>

Developed mainly for linux, but should work for Windows without any major changes. 
