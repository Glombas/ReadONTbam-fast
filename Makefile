# Build readontbam_fast. Needs a C++17 compiler and htslib (>= 1.10) headers + library.
#   make                                  # htslib in a default system location
#   make HTSLIB_PREFIX=$CONDA_PREFIX      # htslib from a conda env (conda install -c bioconda htslib)
#   make HTSLIB_PREFIX=/path/to/htslib    # htslib built from source (its make install prefix)
#   make install PREFIX=$HOME/.local      # copy the binary to PREFIX/bin
CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wno-unused-parameter
HTSLIB_PREFIX ?= $(CONDA_PREFIX)
PREFIX   ?= /usr/local

ifneq ($(strip $(HTSLIB_PREFIX)),)
  CPPFLAGS += -I$(HTSLIB_PREFIX)/include
  LDFLAGS  += -L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
endif
LDLIBS += -lhts -lpthread

BIN = bin/readontbam_fast

all: $(BIN)

$(BIN): src/readontbam_fast.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

install: $(BIN)
	install -d $(PREFIX)/bin
	install -m 755 $(BIN) $(PREFIX)/bin/

clean:
	rm -f $(BIN)

.PHONY: all install clean
