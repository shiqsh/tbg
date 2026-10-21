#!/bin/zsh
# Run from anywhere. Pass program options after the script name.
set -eu
set -o pipefail
task_root=${0:A:h:h}
cd "$task_root/fortran"

# Parallelize k points, not the nested BLAS calls or nested OpenMP regions.
export OMP_NUM_THREADS=${TBG_THREADS:-10}
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_WAIT_POLICY=PASSIVE
export VECLIB_MAXIMUM_THREADS=1

print "OpenMP threads=$OMP_NUM_THREADS; nested OpenMP=1; vecLib threads=1"
# -i prevents idle system sleep while allowing the display to turn off.
# frt is a shell function, so invoke fpm directly here.
exec /usr/bin/caffeinate -i fpm run --target tbg_strain_chi \
    --runner '/usr/bin/time -l' -- "$@"
