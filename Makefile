BOOST = /usr/include/boost/
SEQAN = /env/cns/src/seqan-src/include/

CLASS_FILE=DnaDictionary.cpp ReadFile.cpp SSRContig.cpp SSRContigList.cpp SSRContigLists.cpp NetEx.cpp GeneModel.cpp
CLASS_FILE2=DnaDictionary.cpp ReadFile.cpp
OBJ = ${CLASS_FILE:.cpp=.o}
OBJ2 = ${CLASS_FILE2:.cpp=.o}
OBJ2 += dust.o
MAIN_FILE1=gmorse.cpp
BIN1=gmorse
MAIN_FILE2=kfir.cpp
BIN2=kfir
MAIN_FILE3=gmove.cpp
BIN3=gmove

COMPILER=g++
CPPFLAGS += -I$(BOOST) -I. -I$(SEQAN)
CPPFLAGS += -O3
CPPFLAGS += -pedantic -W -Wall -Wno-long-long -Wno-variadic-macros
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -g
LDFLAGS = -L. -lgzstream -lz
AR = ar cr

COMPILER_BIS=gcc
CFLAGS += -I.
CFLAGS += -O3
CFLAGS += -c

all: $(OBJ)
	${COMPILER} ${CPPFLAGS} -o $(BIN1) $(MAIN_FILE1) $(OBJ)

kfir: libgzstream.a $(OBJ2)
	${COMPILER} ${CPPFLAGS} -o $(BIN2) $(MAIN_FILE2) $(OBJ2) ${LDFLAGS}

gmove: $(OBJ)
	${COMPILER} ${CPPFLAGS} -o $(BIN3) $(MAIN_FILE3) $(OBJ)


gzstream.o:
	${COMPILER} ${CPPFLAGS} -c -o gzstream.o gzstream.cpp

libgzstream.a: gzstream.o
	${AR} libgzstream.a gzstream.o

dust.o:
	${COMPILER_BIS} ${CFLAGS} -o dust.o dust.c

clean:	
	rm -rf *.o a.out 
