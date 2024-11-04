#
# Copy path to external libs here (or comment and set environment variables otherwise)
#
# module load ... # load some modules here if needed on your system
#
# "lib" or "include" end dummy path if required at end of path
export MPI_LIB="/path/to/mpi/lib"
export BOOSTPATH="/path/to/boost/include"
export HDF5_INC="/path/to/hdf5/include"
export HDF5_LIB="/path/to/hdf5/lib"
export SUNDIALS_INC="/path/to/sundials/include"
export SUNDIALS_LIB="/path/to/sundials/lib"
export EIGEN_INC="/path/to/eigen"
export TENSORFLOW_CAPI_PATH="/path/to/tensorflow"
export OPENCV_PATH="/path/to/opencv/folder"
#
# Do not touch from this line
#
script_path="$(dirname -- "${BASH_SOURCE[0]}")"            # relative
script_path="$(cd -- "$MY_PATH" && pwd)"    # absolutized and normalized
echo $script_path
export ORCH_BASE="$script_path/"
export ORCH_MASTER="${ORCH_BASE}"
export ORCH="${ORCH_MASTER}ORCh"
export ORCH_TEST="${ORCH_MASTER}Tests"
export GTCOMB_CT_HOME="${ORCH_BASE}Cantera"
export GTCOMB_CT_HOSTTYPE="${GTCOMB_CT_HOME}/user/linux"
export GTCOMB_CT_DATA="${GTCOMB_CT_HOME}/data"
#
sed -i "/^prefix =/c\prefix = '"${ORCH_BASE}"Cantera_lib'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^python_prefix =/c\python_prefix = '"${ORCH_BASE}"dummy_python'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^python3_prefix =/c\python3_prefix = '"${ORCH_BASE}"dummy_python'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^sundials_include =/c\sundials_include = '"${SUNDIALS_INC}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^sundials_libdir =/c\sundials_libdir = '"${SUNDIALS_LIB}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^boost_inc_dir =/c\boost_inc_dir = '"${BOOSTPATH}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^extra_inc_dirs =/c\extra_inc_dirs = '"${EIGEN_INC}"'" ${GTCOMB_CT_HOME}/cantera.conf
#
export CT_INC="${ORCH_BASE}Cantera_lib/include"
export CT_LIB="${ORCH_BASE}Cantera_lib/lib"
#
export LD_LIBRARY_PATH="${ORCH_BASE}:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${SUNDIALS_LIB}:${LD_LIBRARY_PATH}"
export LIBRARY_PATH="${LIBRARY_PATH}:${OPENCV_PATH}/lib"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${OPENCV_PATH}/lib"
export LIBRARY_PATH="${LIBRARY_PATH}:${TENSORFLOW_CAPI_PATH}/lib"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TENSORFLOW_CAPI_PATH}/lib"
