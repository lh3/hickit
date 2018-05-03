CFLAGS=		-g -Wall -O2 -Wc++-compat -ffast-math
CPPFLAGS=
INCLUDES=
OBJS=		map.o pair.o tad.o neighbor.o phase.o bin.o fdg.o image.o
PROG=		hickit
LIBS=		-lm -lz
ASAN_FLAG=

ifneq ($(asan),)
	ASAN_FLAG = -fsanitize=address
endif

.PHONY:all clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

hickit:$(OBJS) main.o
		$(CC) -o $@ $^ $(ASAN_FLAG) $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *.a *.dSYM hickit.aux hickit.log hickit.pdf

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

bin.o: hkpriv.h hickit.h krng.h khash.h ksort.h
fdg.o: hkpriv.h hickit.h krng.h ksort.h kavl.h
image.o: hkpriv.h hickit.h krng.h ksort.h stb_image_write.h
main.o: hickit.h krng.h
map.o: hickit.h krng.h hkpriv.h khash.h kseq.h
neighbor.o: hkpriv.h hickit.h krng.h ksort.h
pair.o: hkpriv.h hickit.h krng.h ksort.h
phase.o: hkpriv.h hickit.h krng.h
tad.o: hkpriv.h hickit.h krng.h klist.h kavl.h
