SRCS=$(wildcard *.c)
HDRS=$(wildcard *.h)

jhash: $(SRCS) $(HDRS)
	$(CC) -Wall -Wextra -g -O2 -o jellyhash $(SRCS)

clean:
	rm -rf jellyhash *.dSYM

.PHONY: clean
