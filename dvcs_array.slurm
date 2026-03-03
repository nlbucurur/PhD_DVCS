#!/bin/tcsh
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --job-name=dvcs_strip
#SBATCH --array=0-2
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --gres=disk:100000
#SBATCH --output=/w/hallb-scshelf2102/clas12/nlbucuru/PhD_DVCS/logs/%x_%A_%a.out
#SBATCH --error=/w/hallb-scshelf2102/clas12/nlbucuru/PhD_DVCS/logs/%x_%A_%a.err

# # --------------------------
# # Make sure "module" exists
# # --------------------------
# if ( -e /etc/profile.d/modules.csh ) then
#   source /etc/profile.d/modules.csh
# else if ( -e /usr/share/Modules/init/tcsh ) then
#   source /usr/share/Modules/init/tcsh
# endif

# Force temp + caches away from /scratch/slurm
set home = /u/home/$USER
setenv HOME "$home"

if ( -d "$home/tmpfs" ) then
  setenv TMPDIR "$home/tmpfs"
else
  setenv TMPDIR "/tmp/$USER"
endif
mkdir -p "$TMPDIR"

setenv XDG_CACHE_HOME "$TMPDIR/.cache"
mkdir -p "$XDG_CACHE_HOME"
module purge
module use /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/modulefiles
module load clas12
module load coatjava/13.5.2

set BASE = /w/hallb-scshelf2102/clas12/nlbucuru/PhD_DVCS
mkdir -p "$BASE/logs"
mkdir -p "$BASE/farm_outputs"

# tcsh arrays are 1-based, but your Slurm array is 0-2
set PERIODS = ( fall2019 spring2019 spring2020 )
@ INDEX = $SLURM_ARRAY_TASK_ID + 1
set PERIOD = $PERIODS[$INDEX]

set CHANNEL = pDVCS

# Work dir: prefer SLURM_TMPDIR, else TMPDIR (tmpfs/auto), else /tmp
set WORK = "$BASE/work/dvcs_${SLURM_JOB_ID}_${PERIOD}"
mkdir -p "$WORK"
cd "$WORK" || exit 2

echo "Node: `hostname`"
if ( $?SLURM_TMPDIR ) then
  echo "SLURM_TMPDIR = $SLURM_TMPDIR"
else
  echo "SLURM_TMPDIR is not set"
endif

if ( $?TMPDIR ) then
  echo "TMPDIR = $TMPDIR"
else
  echo "TMPDIR is not set"
endif

echo "WORK = $WORK"

# Copy macros/headers needed by includes
cp "$BASE/analysis_DVCS_preskimmed_fiducials_NID.C" .
cp "$BASE/analysis_DVCS_preskimmed_fiducials_NID.h" .
cp "$BASE/run_DVCS.C" .

# Copy any preskimmed helper macros if present (ignore if none)
set nonomatch
cp "$BASE"/*preskimmed*.C . >& /dev/null
unset nonomatch

# Run ROOT non-interactively
cat > run_job.C << EOF
void run_job() {
  gROOT->ProcessLine(".L analysis_DVCS_preskimmed_fiducials_NID.C");
  gROOT->ProcessLine(".L run_DVCS.C");
  gROOT->ProcessLine("run_DVCS(\"$CHANNEL\",\"$PERIOD\");");
}
run_job();
EOF

clas12root -l -b -q run_job.C

if ( $status != 0 ) then
  echo "ERROR: clas12root failed for PERIOD=$PERIOD CHANNEL=$CHANNEL"
  exit 1
endif

# Move outputs back (do not crash if pattern does not match)
set nonomatch
set OUTPAT = stripped_data_${PERIOD}_${CHANNEL}*.root
ls $OUTPAT >& /dev/null
if ( $status == 0 ) then
  mkdir -p "$BASE/farm_outputs/$PERIOD"
  mv -v $OUTPAT "$BASE/farm_outputs/$PERIOD/"
  echo "Wrote: $BASE/farm_outputs/$PERIOD/$OUTPAT"
else
  echo "WARNING: No output matched: $OUTPAT"
endif
unset nonomatch