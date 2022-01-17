#!/usr/bin/env bash
#SBATCH --job-name=PLS_9_13_1
#SBATCH --gres=gpu:1
#SBATCH --qos=qos_gpu-t4
#SBATCH --cpus-per-task=5
#SBATCH --output=./PLS_9_13_1.out
#SBATCH --error=./PLS_9_13_1.err
#SBATCH --time=1500:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
srun python3 Find_PL_Spheres_linear_alg.py