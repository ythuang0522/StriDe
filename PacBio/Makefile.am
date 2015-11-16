noinst_LIBRARIES = libpacbio.a

libpacbio_a_CPPFLAGS = \
	-I$(top_srcdir)/Algorithm \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/Concurrency \
	-I$(top_srcdir)/FMIndexWalk \
	-I$(top_srcdir)/Thirdparty/rollinghash 

libpacbio_a_SOURCES = \
	RollingNode.h RollingNode.cpp \
	RollingPBSelfCTree.h RollingPBSelfCTree.cpp \
	SAIPBSelfCTree.h SAIPBSelfCTree.cpp \
	SAIPBHybridCTree.h SAIPBHybridCTree.cpp \
	PacBioCorrectionProcess.h PacBioCorrectionProcess.cpp
