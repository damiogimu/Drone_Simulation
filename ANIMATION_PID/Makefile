CC=gcc

PLT_FILE1=line_anime.c

P_OBJ1=line_anime.o

PLOT1=lpath
PLOT2=wl

$(PLOT1): $(P_OBJ1)
	$(CC) $(CFLAGS) -o $(PLOT1) $(PLT_FILE1)
	./$(PLOT1)
	rm -f $(PLOT1)

$(PLOT2): $(P_OBJ1)
	$(CC) $(CFLAGS) -D LOAD_F=1 -o $(PLOT2) $(PLT_FILE1)
	./$(PLOT2)
	rm -f $(PLOT2)
