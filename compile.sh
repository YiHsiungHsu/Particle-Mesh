#usage: sh compile.sh filename.cpp

# file name
IN=${1}

# your compiler
CC=g++

# openmp flag of your compiler 

#OPENFLAG=-openmp
#OPENFLAG=-qopenmp
FFTW_PATH=/work1/jared/fftw

# set the output name
OUT=$(basename ${IN} .cpp).out

${CC} -o ${OUT} ${IN} -L${FFTW_PATH}/lib -I${FFTW_PATH}/include  -lfftw3


#g++ -o play play.cpp -L/work1/jared/fftw/lib -I/work1/jared/fftw/include -lfftw3
