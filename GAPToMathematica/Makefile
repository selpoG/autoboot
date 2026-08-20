PROJ := GAPToMathematica
DOTNET ?= dotnet
CONFIGURATION ?= Release
SUFFIXES := .deps.json .dll .pdb .runtimeconfig.json
TARGET := $(addprefix bin/$(PROJ),$(SUFFIXES))
SCRIPTDIR := Scripts
EXEFILE := bin/$(PROJ).dll

SRC := \
	GAPToMathematica.csproj \
	../global.json \
	Makefile \
	*.cs \
	Group/*.cs \
	Parser/*.cs

SCRIPTS := $(wildcard $(SCRIPTDIR)/*.sh)
SCRIPTSTARGET := $(addprefix bin/,$(notdir $(SCRIPTS)))

.PHONY: all run clean

all: $(TARGET) $(SCRIPTSTARGET) ;

$(TARGET) &: $(SRC)
	$(DOTNET) publish -c $(CONFIGURATION) -o bin

bin/%.sh: $(SCRIPTDIR)/%.sh
	cp $< $@
	chmod +x $@

run: all
	$(DOTNET) $(EXEFILE)

clean:
	rm -rf bin obj
