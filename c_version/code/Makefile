CC = gcc
BASE_FLAGS = -Wall -Werror -std=c17 -D _DEFAULT_SOURCE
RELEASE_FLAGS = -O3 $(BASE_FLAGS)
DEBUG_FLAGS = -ggdb $(BASE_FLAGS)
OUTPUT = bin
EXEC = GeneZip
OBJS = $(OUTPUT)/classifier.o $(OUTPUT)/File.o $(OUTPUT)/getline.o $(OUTPUT)/lz78.o $(OUTPUT)/main.o $(OUTPUT)/cli.o
LINK = -lm -pthread

release: COMP_FLAGS = $(RELEASE_FLAGS)
release: clean
release: build

debug: COMP_FLAGS = $(DEBUG_FLAGS)
debug: clean
debug: build

build: clean
build: $(EXEC)

clean:
	rm -f $(OBJS) $(OUTPUT)/$(EXEC)

run: $(EXEC)
	./$(OUTPUT)/$(EXEC)

$(EXEC): $(OBJS)
	$(CC) $^ -o $(OUTPUT)/$(EXEC) $(LINK)

$(OUTPUT)/%.o: %.c
	$(CC) $(COMP_FLAGS) -o $@ -c $?
