LIBS=`root-config --libs`
GLIBS=`root-config --glibs`
CFLAGS=`root-config --cflags`
#CC=/bfactory/package/gcc-2.95.2/bin/g++
EXTRALIBS= "$(ROOTSYS)libMathMore.so"
CC=g++

# set compiler options: 
#   -g  = debugging
#   -O# = optimisation
COPT=-g -c

all: generate_data.exe get_features.exe find_mpv_NN.exe

generate_data.exe:
	$(CC) -o $@ $^ $(CFLAGS) -I. -I$ROOTSYS/lib/ generate_data.C $(GLIBS) $(LIBS) $(CFLAGS) -lpthread -g -ggdb

get_features.exe:
	$(CC) -o $@ $^ $(CFLAGS) -I. -I$ROOTSYS/lib/ get_features.C $(GLIBS) $(LIBS) $(CFLAGS) -lpthread -g -ggdb

find_mpv_NN.exe:
	$(CC) -o $@ $^ $(CFLAGS) -I. -I$ROOTSYS/lib/ find_mpv_NN.C $(GLIBS) $(LIBS) $(CFLAGS) -lpthread -g -ggdb
