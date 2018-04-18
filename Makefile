CFLAGS=		-g -Wall -O2 -Wc++-compat #-fsanitize=address
CPPFLAGS=
INCLUDES=
OBJS=		map.o pair.o tad.o neighbor.o graph.o
PROG=		hickit
LIBS=		-lm -lz

.PHONY:all clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

hickit:$(OBJS) main.o
		$(CC) -o $@ $^ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *.a *.dSYM hickit.aux hickit.log hickit.pdf

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

graph.o: hkpriv.h hickit.h ksort.h
main.o: hickit.h
map.o: hickit.h hkpriv.h khash.h kseq.h
neighbor.o: hkpriv.h hickit.h ksort.h
pair.o: hkpriv.h hickit.h ksort.h
tad.o: hkpriv.h hickit.h klist.h kavl.h
