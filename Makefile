ifdef DEBUG
	OPT=-g -O0 -ggdb
else
	OPT=-O2
endif

all: speedtest

speedtest: speedtest.c jellyhash.h
	$(CC) -Wall -Wextra $(OPT) -o speedtest speedtest.c

clean:
	rm -rf speedtest *.dSYM

.PHONY: all clean
