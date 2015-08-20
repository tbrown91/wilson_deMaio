#!/bin/bash
for i in `seq 1 10`;
do
	./simmlst -N 100 -R 100000 -B 1000000 -s $i
done
echo "Recomb = 10%"
