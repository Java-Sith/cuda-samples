SRC := $(wildcard *.cpp)
EXE := $(SRC:.cpp=)
CFLAGS := -Wall -std=c++17 -O3 -g 
LDFLAGS := -pthread

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC)
	$(CXX) -o $@ $< $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) $(EXE)
