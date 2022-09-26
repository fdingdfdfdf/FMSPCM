# .bash_profile
# Get the aliases and functions
#if [ -f ~/.bashrc ]; then
#	. ~/.bashrc
#fi
# User specific environment and startup programs
#PATH=$PATH:$HOME/bin
#export PATH
alias l='ls -lrt'
alias sq='squeue --user fdingdfdfdf -r'
alias srun='srun -p test --pty --mem 4000 -t 0-08:00  /bin/bash'
alias sa='sacct -j'
alias si='sinfo -o "%15P %5D %4c %7m %10l %20C %N" -e |grep -v _requeue'

#export CLICOLOR=1
#export LSCOLORS=ExFxBxDxCxegedabagacad
#export TERM="xterm-color"

#export exp='/home/fding/fms_web/exp'
#export tmp='/home/fding/scratch-midway/fms_tmp'
#export fv='/n/home01/fdingdfdfdf/cubed_sphere_public_release/exp'
export fv='/n/home01/fdingdfdfdf/cubed_sphere_lbl/exp'
export scratch='/n/regal/wordsworth_lab/fding'
export tmp='/n/regal/wordsworth_lab/fding/fms_tmp'

#export PYTHONPATH=/home/fding/test/lib/python2.7/site-packages/

#module load intel/13.0.079-fasrc01 netcdf/3.6.3-fasrc01
#module load openmpi/1.8.1-fasrc01

#module load centos6/0.0.1-fasrc01

# old FV FMS
#module load intel/13.0.079-fasrc01 netcdf/3.6.3-fasrc01
#module load openmpi/1.8.1-fasrc01

# FMS FVM 
#module load intel/15.0.0-fasrc01
#module load openmpi/1.8.3-fasrc02
#module load netcdf/4.3.2-fasrc04
#module load netcdf-fortran/4.4.0-fasrc02
#module load ncview
#module load nco
#alias ncview='ncview -no_auto_overlay'

#set NETCDF_HOME = /n/sw/fasrcsw/apps/MPI/intel/15.0.0-fasrc01/openmpi/1.8.3-fasrc02/netcdf/4.3.2-fasrc04
#set MPI_HOME = /n/sw/fasrcsw/apps/Comp/intel/15.0.0-fasrc01/openmpi/1.8.3-fasrc02
#PATH=$PATH:$HOME/bin
#export PATH

# Hanieh's options for FMS FVM
#module load intel/13.0.079-fasrc01
#module load openmpi/1.8.1-fasrc01
#module load netcdf/3.6.3-fasrc01
#module load netcdf-fortran/4.4.0-fasrc01

#module load ncview
#module load nco

#alias ncview='ncview -no_auto_overlay'

#module load Anaconda/5.0.1-fasrc02
#conda create -n jupyter_2_7 python=2.7 jupyter 
source activate jupyter_2_7
#export myport=12345
for myport in {6818..11845}; do ! nc -z localhost ${myport} && break; done
echo "ssh -NL $myport:$(hostname):$myport $USER@login.rc.fas.harvard.edu"
jupyter-notebook --no-browser --port=$myport --ip='0.0.0.0'
