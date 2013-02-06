## Process this file with automake to produce Makefile.in

# Detect configuration
if WITH_INST_CTA
INST_CTA = cta
endif
if WITH_INST_LAT
INST_LAT = lat
endif
if WITH_INST_MWL
INST_MWL = mwl
endif
if WITH_INST_COM
INST_COM = com
endif
if WITH_INST_SPI
INST_SPI = spi
endif

# Set sub directories
SUBDIRS = $(INST_MWL) $(INST_CTA) $(INST_LAT) $(INST_COM) $(INST_SPI)

