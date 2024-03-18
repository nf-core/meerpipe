#!/bin/bash

# Check if at least one argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: tempo2_wrapper.sh <arg1> [arg2] [arg3]"
    exit 1
fi

# Access and use the input arguments
archive=$1
ephemeris=$2
select_logic=$3


# Set default value for the third argument if not provided
if [ $# -lt 3 ]; then
    # No third argument provided so don't use -select option
    select_command=""
else
    echo $select_logic > temp_toa_logic.select
    select_command="-select temp_toa_logic.select"
fi

echo $select_logic > temp_toa_logic.select

echo "Generating residuals"
echo "----------------------------------"
# general2 is a tempo2 pluggin: https://github.com/mattpitkin/tempo2/blob/master/plugin/general2_plug.C
# It allows the output of the residuals to be formatted as:
# bat: barycentric arrival time (MJD)
# post: postfit residual (s)
# err: postfit residual error (ns)
# freq: frequency (MHz)
# post_phase: postfit in phase
# flags: flags from the toa files

# Remove lines with 0 snr
sed -i '/ -snr 0.0 /d' ${archive%%.tim}.tim
tempo2 -nofit -output general2 -s "{bat} {post} {err} {freq} {post_phase} {flags}\n" -outfile ${archive}.residual -set START 40000 -set FINISH 99999  -nobs 1000000 -npsr 1 ${select_command} -f ${ephemeris} ${archive%%.tim}.tim  && returncode=$? || returncode=$?
if [[ ${returncode} -ne 134 && ${returncode} -ne 137 && ${returncode} -ne 0 ]]; then
    echo "Errorcode: ${returncode}. Tempo error other than lack of high S/N data error."
    exit $returncode
elif [[ ${returncode} == 134 || ${returncode} == 137 ]]; then
    echo "Errorcode: ${returncode}. No input data due to the logic ${select_logic}"
fi

rm temp_toa_logic.select
