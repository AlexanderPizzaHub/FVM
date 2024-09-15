CC          = g++ -std=c++17

#directory for executables
BIN = bin
PRGM = $(BIN)/DD1D
DEPS = $(shell find ./src -name "*.hpp")
SRC = $(shell find ./src -name "*.cpp")
OBJ = $(SRC:%.c=%.o) 

#--------------------------------------------------
#compiling
#--------------------------------------------------
all: checkdir $(PRGM)

#mkdir
checkdir:
	@mkdir -p $(BIN)

#build executables
%.o: %.c $(DEPS)
	$(CC) -c $< -o $@

$(PRGM): $(OBJ)
	$(CC) $^ -o $@

#clean
.PHONY: clean
clean:
	rm -f $(BIN)/*
	rm ./result/phi*