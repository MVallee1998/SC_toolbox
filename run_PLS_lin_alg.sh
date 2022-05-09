#!/usr/bin/env bash
#SBATCH --job-name=CSPLS_8_12
#SBATCH --gres=gpu:1
#SBATCH --qos=qos_gpu-t4
#SBATCH --cpus-per-task=5
#SBATCH --output=./CSPLS_8_12.out
#SBATCH --error=./CSPLS_8_12.err
#SBATCH --time=1500:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
srun python3 Find_PLS_lin_alg_8_12.py