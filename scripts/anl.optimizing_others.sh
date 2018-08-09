export n_nodes=1
export job_size=$(($COBALT_JOBSIZE / $n_nodes))
export n_mpi_ranks_per_node=1
export n_mpi_ranks=$(($n_nodes * $n_mpi_ranks_per_node))

echo $COBALT_JOBSIZE

export n_hyperthreads_per_core=1

export I_MPI_PIN=1
export I_MPI_DEBUG=4
export KMP_BLOCKTIME=1
export KMP_AFFINITY=scatter
export ARMCI_OPENIB_DEVICE=mlx4_0 
export n_openmp_threads_per_rank=64
#export n_hyperthreads_skipped_between_ranks=64
export MPICH_GNI_MAX_VSHORT_MSG_SIZE=10000
export MPICH_GNI_MAX_EAGER_MSG_SIZE=98304
export MPICH_GNI_NUM_BUFS=300
export MPICH_GNI_NDREG_MAXSIZE=16777216
export MPICH_GNI_MBOX_PLACEMENT=nic
export COMEX_MAX_NB_OUTSTANDING=6 

unset OMP_NESTED
unset MKL_DYNAMIC
unset KMP_PLACE_THREADS
unset OMP_DYNAMIC
unset MKL_NUM_THREADS
export MKL_DYNAMIC=false
export OMP_DYNAMIC=false
export OMP_NESTED=true


export OMP_PLACES="cores"
export OMP_PROC_BIND=spread,close
export I_MPI_PIN_DOMAIN=auto


# see aprun --help for more options
#--env FOO=$FOO --env BAR=$BAR

for i in `seq 1 $job_size`;
    do
        echo $i
        aprun -n $n_mpi_ranks -N $n_mpi_ranks_per_node \
            --env OMP_NUM_THREADS=$n_openmp_threads_per_rank -cc depth \
            --env I_MPI_PIN=1 \
            --env KMP_BLOCKTIME=1 \
            --env KMP_AFFINITY=scatter \
            --env ARMCI_OPENIB_DEVICE=mlx4_0 \
            --env MPICH_GNI_MAX_VSHORT_MSG_SIZE=10000  \
            --env MPICH_GNI_MAX_EAGER_MSG_SIZE=98304  \
            --env MPICH_GNI_NUM_BUFS=300 \
            --env MPICH_GNI_NDREG_MAXSIZE=16777216  \
            --env MPICH_GNI_MBOX_PLACEMENT=nic \
            --env COMEX_MAX_NB_OUTSTANDING=6 \
            python ~/Code/conformer/scripts/optimizing_others.py $i > \
                /projects/CPOX/northeastern_comocheng/conformers/optimizing_others/others.$i.log &
        
    done 