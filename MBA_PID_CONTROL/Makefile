CC=gcc
CFLAGS=-Wall

HEADER=my_header.h
OBJS=sim_mba_pid_control.o utils.o
P_OBJ1=DRAW_FILES/draw_error.o
P_OBJ_ACC=DRAW_FILES/draw_acc.o
P_OBJ_VELO=DRAW_FILES/draw_velo.o
P_OBJ_POS=DRAW_FILES/draw_pos.o
P_CFILE3=DRAW_FILES/draw_data.c
P_OBJ3=DRAW_FILES/draw_data.o

CALC=calc
PLOT1=lerror
PLOT_ACC=lacc
PLOT_VELO=lvelo
PLOT_POS=lpos
PLOT_ATTITUDE=latti
PLOT_PAYLOAD=lpayload

$(CALC): $(OBJS)
	$(CC) -o $@ $^ -lm	
	./$@
	rm $@

$(PLOT1): $(P_OBJ1)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT_ACC): $(P_OBJ_ACC)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT_VELO): $(P_OBJ_VELO)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT_POS): $(P_OBJ_POS)
	$(CC) -o $@ $^
	./$@
	rm -f $@

$(PLOT_ATTITUDE): $(P_CFILE3)
	$(CC) -D ATTI_F=1 -o $@ $^
	./$@
	rm -f $@

$(PLOT_PAYLOAD): $(P_OBJ3)
	$(CC) -o $@ $^
	./$@
	rm -f $@

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $< 
