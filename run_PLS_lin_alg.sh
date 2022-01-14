#!/usr/bin/env bash
#SBATCH --job-name=MyFirstJob
#SBATCH --gres=gpu:2
#SBATCH --qos=qos_gpu-t4
#SBATCH --cpus-per-task=5
#SBATCH --output=./MyFirstJob.out
#SBATCH --error=./MyFirstJob.err
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
srun python3 Find_PL_Spheres_linear_alg.py