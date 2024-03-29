##########################################
# SET THESE 6 PATHS CORRECTLY TO COMPILE #
##########################################
BOOST_INC=
BOOST_LIB=
HTSLD_INC=
HTSLD_LIB=
#########################################################
# EXAMPLES:                                             #
# BOOST_INC=/usr/include/                               #
# BOOST_LIB=/usr/lib/x86_64-linux-gnu/                  #
# HTSLD_INC=$(HOME)/Tools/htslib-1.9                    #
# HTSLD_LIB=$(HOME)/Tools/htslib-1.9                    #
#########################################################


define n


endef

#INSTALL LOCATIONS
#CHANGE prefix TO INSTALL LOCALLY
prefix=/Users/nikolaoslykoskoufis/Documents/PROJECTS/fenrichcpp
exec_prefix=$(prefix)
bindir=$(exec_prefix)/bin
datarootdir=$(prefix)/share
MKDIR_P=mkdir -p
INSTALL=install -p
INSTALL_EXE=$(INSTALL) -m 755
INSTALL_DIR=$(MKDIR_P) -m 755
INSTALL_PROGRAM=$(INSTALL)
INSTALL_SCRIPT=$(INSTALL_EXE)


#COMPILER MODE C++11 
CXX=g++ -std=c++0x

#COMPILER FLAGS 
CXXFLAG_REL=-O3
CXXFLAG_DBG=-g
CXXFLAG_WRN=-Wall -Wextra -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated -Wno-unused-parameter


#BASE LIBRARIES
LIB_FLAGS=-lz -lgsl -lbz2 -llzma -lgslcblas -lm -lpthread -lcurl 

#FILE LISTS 
BFILE=bin/fenrich
HFILE=$(shell find src -name *.h)
TFILE=$(shell find lib -name *.h)
OFILE=$(shell for file in `find src -name *.cpp | LC_ALL=C sort`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp | LC_ALL=C sort`; do echo $$(dirname $$file); done)



#STATICLY LINKED LIBS
LIB_FILES=$(HTSLD_LIB)/libhts.a $(BOOST_LIB)/libboost_iostreams.a $(BOOST_LIB)/libboost_program_options.a

#INCLUDE DIRS
IFLAG=-I$(HTSLD_INC) -I$(BOOST_INC) -isystem lib

# FOR MAC ONLY 
MZ=/usr/local/opt/zlib/lib/libz.a
MCBLAS=/usr/local/lib/libgslcblas.a
MGSL=/usr/local/lib/libgsl.a
MBZ2=/usr/local/opt/bzip2/lib/libbz2.a
MLZMA=/usr/local/lib/liblzma.a
MCURL=/usr/local/opt/curl/lib/libcurl.a

#MAC SPECIFIC STUFF
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
    CXXFLAG_REL+= -fvisibility=hidden -fvisibility-inlines-hidden
    CXXFLAG_DBG+= -fvisibility=hidden -fvisibility-inlines-hidden
    autocompdir = /usr/local/etc/bash_completion.d
 ifeq ($(MAKECMDGOALS),static)
    ifeq ("$(wildcard  $(MZ))","")
    $(error Cannot find $(MZ)! $nTry:	brew install zlib$nOr edit the makefile with the correct    location of libz.a)
   endif
   ifeq ("$(wildcard  $(MCBLAS))","")
    $(error Cannot find $(MCBLAS)! $nTry:	brew install gsl$nOr edit the makefile with the correct location of libgslcblas.a)
   endif
   ifeq ("$(wildcard  $(MGSL))","")
    $(error Cannot find $(MGSL)! $nTry:	brew install gsl$nOr edit the makefile with the correct location of libgsl.a)
   endif
   ifeq ("$(wildcard  $(MBZ2))","")
    $(error Cannot find $(MBZ2)! $nTry:	brew install bzip2$nOr edit the makefile with the correct location of libbz2.a)
   endif
   ifeq ("$(wildcard  $(MLZMA))","")
    $(error Cannot find $(MLZMA)! $nTry:	brew install xz$nOr edit the makefile with the correct location of liblzma.a)
   endif
   ifeq ("$(wildcard  $(MCURL))","")
    $(error Cannot find $(MCURL)! $nTry:	brew install curl$nOr edit the makefile with the correct location of libcurl.a)
   endif           
 endif
endif

#CHECK IF WE MADE BEFORE INSTALLING
ifeq ($(MAKECMDGOALS),install)
    ifeq ("$(wildcard $(BFILE))","")
    $(error Cannot find $(BFILE)! $nPlease make before you make install)
    endif
endif

#STATIC VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
static: CXXFLAG=$(CXXFLAG_REL) $(CXXFLAG_WRN)
static: LDFLAG=$(CXXFLAG_REL)
static: $(BFILE)


#STATIC VERSION (SET UP THE VARIABLES IN THE BEGINING OF THE MAKEFILE)
personal: BOOST_INC=/Users/nikolaoslykoskoufis/Documents/Programming/Tools/boost_1_74_0/
personal: BOOST_LIB=/Users/nikolaoslykoskoufis/Documents/Programming/Tools/boost_1_74_0/stage/lib/
personal: HTSLD_INC=/Users/nikolaoslykoskoufis/Documents/Programming/Tools/htslib-1.11/
personal: HTSLD_LIB=/Users/nikolaoslykoskoufis/Documents/Programming/Tools/htslib-1.11/
personal: static

baobab: BOOST_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/
baobab: BOOST_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/boost_1_71_0/stage/lib/
baobab: HTSLD_INC=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab: HTSLD_LIB=/srv/beegfs/scratch/groups/funpopgen/Tools/htslib-1.9/
baobab: static


install:
	$(INSTALL_DIR) $(DESTDIR)$(bindir) 
	$(INSTALL_PROGRAM) $(BFILE)


#COMPILATION RULES 
$(BFILE):$(OFILE)
	$(CXX) $^ $(LIB_FILES) -o $@ $(LIB_FLAGS) $(LDFLAG) $(IFLAG)

obj/fenrich.o: src/fenrich.cpp $(HFILE) $(TFILE) $(CFILE)
	$(CXX) -o $@ -c $< ${CXXFLAG} $(IFLAG)

obj/null_%.o: null_%.cpp null_data.h $(TFILE)
	$(CXX) -o $@ -c $< ${CXXFLAG} $(IFLAG)

obj/analysis_%.o: analysis_%.cpp analysis_data.h  $(TFILE)
	$(CXX) -o $@ -c $< ${CXXFLAG} $(IFLAG)


clean:
	rm -f obj/*.o $(BFILE)

clean-null:
	rm -f obj/*.o $(BFILE)

clean-analysis:
	rm -f obj/*.o $(BFILE)
