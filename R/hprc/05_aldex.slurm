#!/bin/bash
##############################################################################
# Submit 05a (memory-lean ALDEx2) on HPRC.
#
# The whole point of 05a is that peak RAM is predictable:
#
#     GB ~= n_transcripts * n_libraries * mc.samples * 8 / 1e9 * 2.5
#
# Set --mem below from that number with ~20% headroom. 05a prints its own
# prediction and refuses to start if it exceeds mem_budget_gb, so a wrong
# guess here fails in seconds rather than 6 hours in.
#
#   207k transcripts x 23 libraries, mc.samples = 128 -> ~122 GB -> --mem=150G
#                                    mc.samples =  64 ->  ~61 GB -> --mem=80G
#
# Cheapest way to cut this is a stricter prevalence filter, not fewer MC
# instances -- see prevalence_threshold_report.csv written by 05a.
##############################################################################
#SBATCH --job-name=aldex-lean
#SBATCH --output=logs/aldex-lean-%j.out
#SBATCH --error=logs/aldex-lean-%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G

set -euo pipefail
mkdir -p logs

module purge
module load GCC/13.2.0 OpenMPI/4.1.6 R/4.4.1    # adjust to the current stack

echo "host: $(hostname)   started: $(date)"

Rscript 05a_aldex_lean_HPRC.R

echo "finished: $(date)"

# Then, in a much smaller job (or interactively -- it needs a few GB at most):
#   Rscript 05b_wrangle_de_HPRC.R
