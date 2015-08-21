#!/bin/bash
for i in `seq 1 100`;
do
	./simmlst -N 100 -R 0 -B 1000000 -s $i >> out_0
done

for i in `seq 1 100`;
do
	./simmlst -N 100 -R 100 -B 1000000 -s $i >> out_0.01
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 100 -B 1000000 -s $i >> out_0.1
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 1000 -B 1000000 -s $i >> out_1
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 2000 -B 1000000 -s $i >> out_2
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 5000 -B 1000000 -s $i >> out_5
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 10000 -B 1000000 -s $i >> out_10
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 15000 -B 1000000 -s $i >> out_15
done

for i in `seq 1 100`;
do
	./simmlst -N 1000 -R 20000 -B 1000000 -s $i >> out_20
done

