# Author SHI Xin <shixin@ihep.ac.cn>  
# Created [2021-03-31 Sun 14:38] 


BIN=./bin
SRC=./src

CC = g++
GCCFLAGS  = -Wall -g 

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))

FLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lHistPainter 


PROG=raser
LIST=$(addprefix $(BIN)/, $(PROG))


all: $(LIST)
	@echo "Build successful."

$(LIST): | $(BIN)

$(BIN): 
	mkdir -p $(BIN)

$(BIN)/raser: $(SRC)/raser.cc
	$(CC) $< $(FLAGS) -o $@

merge: 
	git remote update && git merge dt-np/main 

clean:
	rm -f $(BIN)/raser 
	