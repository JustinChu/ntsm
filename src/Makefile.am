bin_PROGRAMS = ntfp

ntfp_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

ntfp_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
ntfp_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-lz \
	-pthread

ntfp_SOURCES = ntfp.cpp \
	Options.h \
	Options.cpp \ 
	Util.h