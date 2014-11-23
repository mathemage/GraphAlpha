CC = gcc
CFLAGS = -g -O1
OBJECTS = main.o GraphAlpha.o

twhree : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o GraphAlpha

%.o : %.c
	ctags *.[ch] -R && $(CC) $(CFLAGS) -c $<

clean:
	rm $(OBJECTS)
