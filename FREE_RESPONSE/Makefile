CC=gcc
CFLAGS=-Wall
HEADER=my_header.h

OBJ=sim_noncontrol.o
P_OBJ=draw.o
CALC=calc
PLOT=plot

$(CALC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@

$(PLOT): $(P_OBJ)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@

%.o: %.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
