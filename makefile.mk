CC	=g++
CFLAGS	=-c -Wall -O4 -ffast-math -ffinite-math-only -I header/
LDFLAGS	=
SOURCES	=./examples/test_CUR_Lowrank.cpp
MATRIX	=-DKERNEL  # use -DEXACTLOWRANK, -DKERNEL
ROWCOL	=-DEQUISPACED	# -DEQUISPACED, -DCHEBSPACED
OBJECTS	=$(SOURCES:.cpp=.o)
EXECUTABLE	=./exec/CUR_Lowrank

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(MATRIX) $(ROWCOL) $< -o $@

clean:
	rm -rf *.out ./examples/*.o ./exec/*

tar:
	tar -zcvf CUR_Lowrank.tar.gz ./makefile.mk ./exec ./header ./examples ./README.md ./LICENSE.md