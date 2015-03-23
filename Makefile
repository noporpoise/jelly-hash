ifdef DEBUG
	OPT=-g -O0 -ggdb
else
	OPT=-O3
endif

all: speedtest

speedtest: speedtest.c jellyhash.h num_format.c num_format.h
	$(CC) -Wall -Wextra $(OPT) -o speedtest speedtest.c num_format.c -lpthread -lm

clean:
	rm -rf speedtest *.dSYM

.PHONY: all clean
