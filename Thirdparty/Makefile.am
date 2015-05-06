AUTOMAKE_OPTIONS = subdir-objects
noinst_LIBRARIES = libthirdparty.a

libthirdparty_a_CPPFLAGS = \
        -I$(top_srcdir)/FMIndexWalk \
        -I$(top_srcdir)/Algorithm \
        -I$(top_srcdir)/Bigraph \
        -I$(top_srcdir)/Util \
        -I$(top_srcdir)/StringGraph \
        -I$(top_srcdir)/SQG \
        -I$(top_srcdir)/SuffixTools \
        -I$(top_srcdir)/GraphDifif \
    -I$(top_srcdir)/Concurrency

libthirdparty_a_SOURCES = \
        gzstream.h gzstream.C \
        stdaln.h stdaln.c \
        overlapper.h overlapper.cpp \
        multiple_alignment.h multiple_alignment.cpp \
        bcr.h bcr.c \
	ropebwt2/mrope.h ropebwt2/mrope.c \
	ropebwt2/rle.h ropebwt2/rle.c \
	ropebwt2/rope.h ropebwt2/rope.c \
	ropebwt2/kseq.h \
	ropebwt2/rld0.h ropebwt2/rld0.c
