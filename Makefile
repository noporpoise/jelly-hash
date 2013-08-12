jellytest: jellytest.c jellyhash.h
	$(CC) -Wall -Wextra -g -O2 -o jellytest jellytest.c

clean:
	rm -rf jellytest *.dSYM

.PHONY: clean
