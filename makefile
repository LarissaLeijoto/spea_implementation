# Directories.
export BINDIR    = $(CURDIR)/bin
export SRCDIR    = $(CURDIR)/src

# Executable name.
export EXEC = spea2

# Toolchain.
export CXX=g++ -g -std=c++11

# Toolchain configuration. -std=c99  -Werror
export CXXFLAGS = -I -ansi -pedantic -Wall -Wextra -Wno-sign-compare 
export CXXFLAGS += -D NDEBUG -O3

# Builds everything.
all:
	@mkdir -p $(BINDIR)
	@cd $(SRCDIR) && make all

# Cleans compilation files.
clean:
	@cd $(SRCDIR) && make clean
