GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
CC = gcc
CFLAGS += -g -std=c99 -I$(GSLINC) -I./tools -I./EOBNRv2HMROM -I./integration -I./LISAsim -I./LLVsim -I./LLVinference

SUBDIRS = tools EOBNRv2HMROM LISAsim LLVsim integration LLVinference
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

.PHONY: all clean subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

EOBNRv2HMROM: tools

LISAsim: tools integration EOBNRv2HMROM

LLVsim: tools integration EOBNRv2HMROM

LLVinference: tools integration EOBNRv2HMROM LLVsim

all: subdirs

clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
