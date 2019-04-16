PROJ := GAPToMathematica
SUFFIXES := .deps.json .dll .pdb .runtimeconfig.json
TARGET := $(addprefix bin/$(PROJ),$(SUFFIXES))
SCRIPTDIR := Scripts
EXEFILE := bin/$(PROJ).dll

SRC := \
	*.cs \
	Group/*.cs \
	Parser/*.cs

SCRIPTS := $(wildcard $(SCRIPTDIR)/*.sh)
SCRIPTSTARGET := $(addprefix bin/,$(notdir $(SCRIPTS)))

.PHONY: all run clean

all: $(TARGET) $(SCRIPTSTARGET) ;

$(TARGET): $(SRC)
	dotnet publish -o bin

bin/%.sh: $(SCRIPTDIR)/%.sh
	cp $< $@

run: all
	mono $(EXEFILE)

clean:
	rm -f $(TARGET)
	rm -f $(SCRIPTSTARGET)
	rm -rf bin/Debug
	rm -rf obj
