CC=gcc
CFLAGS=-Wall

HEADER=my_header.h
OBJS=sim_pid_control.o utils.o
P_OBJ1=DRAW_FILES/draw_error.o
P_OBJ2=DRAW_FILES/draw_acc.o
P_CFILE3=DRAW_FILES/draw_data.c
P_OBJ3=DRAW_FILES/draw_data.o

CALC=calc
PLOT1=lerror
PLOT2=lacc
PLOT_ATTITUDE=latti
PLOT_PAYLOAD=lpayload

$(CALC):$(OBJS)
	$(CC) -o $@ $^
	./$@
	rm $@

$(PLOT1):$(P_OBJ1)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT2):$(P_OBJ2)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT_ATTITUDE):$(P_CFILE3)
	$(CC) -D ATTI_F=1 -o $@ $^
	./$@
	rm -f $@

$(PLOT_PAYLOAD):$(P_OBJ3)
	$(CC) -o $@ $^
	./$@
	rm -f $@

clean:
	rm -f $(OBJS)

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
