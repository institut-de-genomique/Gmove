BOOST = /usr/include/boost/
#SEQAN = /env/cns/src/seqan-src/include/

DIR_BIN_PROD=$(DIR)

CLASS_FILE=GffRecord.cpp GffRecordList.cpp SSRContig.cpp SSRContigList.cpp SSRContigLists.cpp NetEx.cpp GeneModel.cpp GeneModelList.cpp
OBJ = ${CLASS_FILE:.cpp=.o}
MAIN_FILE=gmove.cpp
BIN=gmove

COMPILER=g++
CPPFLAGS += -I$(BOOST) -I. #-I$(SEQAN)
CPPFLAGS += -O3
CPPFLAGS += -W -Wall -Wno-long-long -Wno-variadic-macros -Wno-deprecated
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -g
AR = ar cr

COMPILER_BIS=gcc
CFLAGS += -I.
CFLAGS += -O3
CFLAGS += -c


gmove: $(OBJ)
	${COMPILER} ${CPPFLAGS} -o $(DIR_BIN_PROD)$(BIN) $(MAIN_FILE) $(OBJ)

clean:	
	rm -rf *.o a.out 