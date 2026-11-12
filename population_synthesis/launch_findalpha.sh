#!/bin/bash

# Create logs directory if it does not exist
if [ ! -d "logs" ]; then
    mkdir -p logs
fi

# Input names
names=(
    B1534+12
    B1913+16
    J0509+3801
    J0653+0443
    J0737-3039A
    J1208-5936
    J1756-2251
    J1757-1854
    J1846-0513
    J1856-0039
    J1906+0746
    J1913+1102
    J1946+2052
)

names1=(
    B1534+12
    B1913+16
    J0509+3801
    J0737-3039A
    J1208-5936
    J1756-2251
    J1757-1854
    J1856-0039
    J1913+1102
    J1946+2052
)

names2=(J1906+0746)

names3=(
    J0653+0443
    J1846-0513
)


# Submit one independent job per name
for name in "${names1[@]}"; do
    sbatch \
        --job-name="find_alpha_${name}_nG" \
        --output="logs/find_alpha_${name}_nG.out" \
        --error="logs/find_alpha_${name}_nG.err" \
        findalpha.sbatch "$name"
done

