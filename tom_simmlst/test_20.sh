#!/bin/bash
for i in `seq 1 1`;
do
	./simmlst -N 100 -R 400000 -B 2000000 -s $i
done
echo "Recomb = 20%, genome = 2MB"
