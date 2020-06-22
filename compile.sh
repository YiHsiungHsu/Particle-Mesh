#usage: sh compile.sh filename.cpp

# filename
IN=${1}

# your compiler
CC=mpic++

# openmp flag of your compiler
OPENFLAG=-fopenmp

#fftw path
FFTW_PATH=/work1/jared/fftw

# set the output name
OUT=$(basename ${IN} .cpp).out

${CC} ${OPENFLAG} ${IN} -o ${OUT} -L${FFTW_PATH}/lib -I${FFTW_PATH}/include  -lfftw3

#g++ -o play play.cpp -L/work1/jared/fftw/lib -I/work1/jared/fftw/include -lfftw3
