PROJ := GAPToMathematica
TARGET := $(PROJ).exe
SCRIPTDIR := Scripts
EXEDIR := bin/Debug

SRC := \
	$(PROJ)/*.cs \
	$(PROJ)/Group/*.cs \
	$(PROJ)/Parser/*.cs \
	$(PROJ)/Properties/*.cs

SCRIPTS := $(wildcard $(PROJ)/$(SCRIPTDIR)/*.sh)
SCRIPTSTARGET := $(notdir $(SCRIPTS))

.PHONY: all run clean

all: $(TARGET) $(SCRIPTSTARGET) ;

$(PROJ)/$(EXEDIR)/$(TARGET): $(SRC)
	xbuild

$(TARGET): $(PROJ)/$(EXEDIR)/$(TARGET)
	cp $(PROJ)/$(EXEDIR)/$(TARGET) $(TARGET)

%.sh: $(PROJ)/$(SCRIPTDIR)/%.sh
	cp $< $@

run: all
	mono $(TARGET)

clean:
	rm -f $(TARGET)
	rm -f $(SCRIPTSTARGET)
	rm -rf $(PROJ)/bin
	rm -rf $(PROJ)/obj
