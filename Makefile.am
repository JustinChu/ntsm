SUBDIRS = src

dist_doc_DATA = README.md

bin_SCRIPTS = ntsmSiteGen

EXTRA_DIST = vendor \
			ntsmSiteGen \
			ntsm-scripts

ntsmScriptsdir = $(prefix)/bin/ntsm-scripts
ntsmScripts_DATA = ntsm-scripts/convertTSVtoPCA.py \
				ntsm-scripts/extractSNPsfromVCF.py \
				ntsm-scripts/filterRepetiveSNP.pl \
				ntsm-scripts/makefile
