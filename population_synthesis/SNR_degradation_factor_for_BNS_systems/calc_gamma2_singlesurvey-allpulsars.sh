#!/bin/bash

for s in 1534 1913 0509 0737A 0737B 1208 1756 1757 1856 1906 1102 1946; do
    cp params_$s.in params.in
    ./gamma2exec >> GPPS.txt
done
