
# THIS MAKEFILE EXISTS ONLY TO BUILD THE PURE C LIBRARY VERSION OF intervaldb.c

#Ignore the Python-dependant stuff in Pygr's default.h
#will tidy this up at some point.
CFLAGS=-DBUILD_C_LIBRARY

INC=intervaldb.h
SRC=intervaldb.c
OBJ=intervaldb.o
OUT=libnclist.so



#not really sure how this works, but it seems to.
all: $(OUT)

$(OUT): $(OBJ)
	$(CC) -shared -o $@ $^

%.o: %.c $(INC)
	$(CC) $(CFLAGS) -fpic -c $< -o $@

clean:
	rm -f $(OBJ) $(OUT)
