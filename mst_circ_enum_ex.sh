#!/bin/bash
python3 /home/morriang/circuit_enumeration/MST_testing/main.py "365" "400"

Outer bash, named: mst_circ_enum.sh 
#!/bin/bash

#SBATCH --job-name=mst_circ_enum

#SBATCH --output=mst_circ_enum.out.%j

#SBATCH --error=mst_circ_enum.err.%j

#SBATCH -n 8

#SBATCH -p math-alderaan-gpu

# singularity exec /storage/singularity/pyscipopt-geopandas.sif /home/hortondr/equitable_facility_location/inner_kpcon_single.sh
