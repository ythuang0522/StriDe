bin_PROGRAMS = stride

stride_CPPFLAGS = \
	-I$(top_srcdir)/Util \
	-I$(top_srcdir)/Bigraph \
	-I$(top_srcdir)/SuffixTools \
	-I$(top_srcdir)/StringGraph \
	-I$(top_srcdir)/Concurrency \
	-I$(top_srcdir)/Algorithm \
	-I$(top_srcdir)/SQG \
	-I$(top_srcdir)/Thirdparty \
	-I$(top_srcdir)/FMIndexWalk

stride_LDADD = \
	$(top_builddir)/StringGraph/libstringgraph.a \
	$(top_builddir)/Concurrency/libconcurrency.a \
	$(top_builddir)/Algorithm/libalgorithm.a \
	$(top_builddir)/SuffixTools/libsuffixtools.a \
	$(top_builddir)/Bigraph/libbigraph.a \
	$(top_builddir)/Util/libutil.a \
	$(top_builddir)/SQG/libsqg.a \
	$(top_builddir)/Thirdparty/libthirdparty.a \
	$(top_builddir)/FMIndexWalk/libfmindexwalk.a

stride_LDFLAGS = -pthread

stride_SOURCES = StriDe.cpp \
		strideall.cpp strideall.h\
              index.cpp index.h \
              overlap.cpp overlap.h \
              assemble.cpp assemble.h \
              correct.cpp correct.h \
              oview.cpp oview.h \
              preprocess.cpp preprocess.h \
              subgraph.cpp subgraph.h \
              filter.cpp filter.h \
              fm-merge.cpp fm-merge.h \
              OverlapCommon.h OverlapCommon.cpp \
		kmerfreq.h kmerfreq.cpp \
		grep.h grep.cpp \
		FMIndexWalk.h FMIndexWalk.cpp \
              SGACommon.h 
