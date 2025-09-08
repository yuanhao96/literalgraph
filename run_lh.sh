#!/bin/bash

#SBATCH --job-name=eval
#SBATCH --nodes=1
#SBATCH --mail-type=END
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=0
#SBATCH --time=2-00:00:00
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-a100
#SBATCH --mem=300g
#SBATCH --output=eval.log
#SBATCH --mail-user=hyhao@umich.edu
#SBATCH --mail-type=BEGIN,END

# conda activate kg
export LD_LIBRARY_PATH=/home/hyhao/.conda/envs/kg/lib:$LD_LIBRARY_PATH
module load cuda/11.8.0 cudnn/11.8-v8.7.0
# module load python/3.9.12
# module load numpy pandas/1.4.2

export OPENAI_API_KEY=**********
python scripts/evaluate_ner.py
