bin_PROGRAMS = ntfp

ntfp_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

ntfp_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
ntfp_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-pthread

ntfp_LDADD = -lz

ntfp_SOURCES = ntfp.cpp