ROOTLIBS = `root-config --glibs`
ROOTCFLAGS = `root-config --cflags`

BIN := ../bin
SRC_DIR := ../src
LDLIBS=$(ROOTLIBS) -lgsl -lgslcblas -L/seaquest/users/cmantill/DarkQuest/lhe/HepMC/lib/ -lHepMC -I/seaquest/users/cmantill/DarkQuest/lhe/HepMC/include/
INCFILES=-I$(SRC_DIR)
CPP=g++
CPPFLAGS=$(LDLIBS) $(ROOTCFLAGS) -g
TARGETS := $(patsubst $(SRC_DIR)/%.cc,$(BIN)/%,$(wildcard $(SRC_DIR)/*.cc))
UTIL_SRC := $(wildcard $(SRC_DIR)/*.cpp)
UTIL_OBJ := $(UTIL_SRC:.cpp=.o)

all: $(TARGETS) $(UTIL_OBJ)
$(BIN)/%: $(SRC_DIR)/%.cc $(UTIL_OBJ)
	$(CPP) $< $(UTIL_OBJ) $(CPPFLAGS) $(INCFILES) -o $@ 
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.hh
	$(CPP) $< -c $(INCFILES) -g -o $@
clean:
	@rm -f $(TARGETS) $(UTIL_OBJ) 

