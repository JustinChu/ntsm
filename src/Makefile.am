bin_PROGRAMS = ntsmCount ntsmEval

ntsmCount_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

ntsmCount_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
ntsmCount_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-pthread

ntsmCount_LDADD = -lz

ntsmCount_SOURCES = ntSeqMatchCount.cpp \
	FingerPrint.hpp \
	Options.cpp Options.h \
	Util.h
	

ntsmEval_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

ntsmEval_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
ntsmEval_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-pthread

ntsmEval_SOURCES = ntSeqMatchEval.cpp \
	CompareCounts.hpp \
	Options.cpp Options.h \
	Util.h
