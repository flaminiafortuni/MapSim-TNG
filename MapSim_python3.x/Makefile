# executable name 
PROG = /YOUR_PATH/MapSim-TNG/TESTMapSim #change ex

MAIN = TESTmain.cpp readSUBFIND.cpp

LDFLAGS +=-Wl,-rpath  -lstdc++ #-DSWAP #-DCOSMOLIB 

##NR_DIR = /marconi/home/userexternal/cgiocoli/lib/NR
COSMOLIB_DIR = /YOUR_PATH/CosmoLib


# gsl, cfitsio, CCfits, fftw
LIBS = -L/usr/local/lib/gsl/  -lgsl -lgslcblas  \
       -L/home/YOUR_USER/.local/lib/python3.6/site-packages/numpy/ \ #if you have conda
       -L/usr/lib/python3.6/ -lpython3.6m \
       -L/usr/local/lib/cfitsio/ \
       -L/usr/local/lib/CCfits/ -lcfitsio -lCCfits \
       -L/usr/local/lib/fftw/ -lfftw3 -lm \
       -L/usr/lib64 -lstdc++
       #-L$(COSMOLIB_DIR)/ -lCosmoLib \
       #-L$(NR_DIR)/ -lNR \

# gsl, cfitsio, CCfits, fftw  
ALLFLAGS = -I/usr/local/include/gsl/ \
	   -I/usr/local/include/gsl/ \
           -I/usr/local/include/cfitsio/ \
           -I/usr/local/include/CCfits/ \
           -I/usr/local/include/fftw/ \
           -I/home/YOUR_USER/.local/lib/python3.6/site-packages/numpy/core/include/numpy/ \
           -I/usr/include/python3.6m/ \
	   -I./
           #-I/usr/include/python2.7/ -lpython2.7 \ #add by me
           #-I/home/flaminia/.local/lib/python2.7/site-packages/numpy/core/include/numpy/ \ #add by me
	   #-I$(NR_DIR)/include \
           #-I$(COSMOLIB_DIR)/include \

# 
DEBUG = -g -O2 

CC = /usr/bin/cc -std=c++11 -lstdc++

#
RM = rm -f -r
#
OBJ = $(SOURCES:.cpp=.o)
#

CFLAGS=-O2  -g -fPIC -DPIC

default: main
main: 
	$(CC) $(CFLAGS) ${ALLFLAGS} $(MAIN) ${LIBS} ${LDFLAGS} -o ${PROG}
clean:
	$(RM) $(PROG) $(OBJ) *~

