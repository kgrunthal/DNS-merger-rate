#!/bin/bash

# PALFA		268s 
# PMSURV	2100s
# AODRIFT	50s	(B1534+12)
# GBNCC		120s
# HTRU-low	340s
# (HTRU-mid	540s)
# PHSURV	265s	
# MMGPS-L	600s
# GPPS		300s	


for t in 268 2100 50 120 340 265 600 300; do
    sed -i "2s/.*/${t}.0D0/" params.in
    ./gamma2exec >> 1856.txt
done
