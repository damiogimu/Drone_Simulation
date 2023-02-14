CC=gcc
CFLAGS=-Wall

HEADER=my_header.h
OBJS=sim_pid_control.o trajectories.o utils.o

P_OBJ_ACC=DRAW_FILES/draw_acc.o
P_OBJ_VELO=DRAW_FILES/draw_velo.o
P_OBJ_POS=DRAW_FILES/draw_pos.o
P_OBJ_ATTITUDE=DRAW_FILES/draw_attitude.o
P_OBJ_LOAD=DRAW_FILES/draw_load.o
3D_GIF_OBJ=DRAW_FILES/3d_anime.o
2D_GIF_OBJ=DRAW_FILES/2d_anime.o

CALC=calc
PLOT_ACC=p_acc
PLOT_VELO=p_velo
PLOT_POS=p_pos
PLOT_ATTITUDE=p_atti
PLOT_PAYLOAD=p_load
PLOT_GIF_3d=3d
PLOT_GIF_2d=2d

$(CALC): $(OBJS)
	$(CC) -o $@ $^ -lm	
	./$@
	rm $@

$(PLOT_ACC): $(P_OBJ_ACC)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_VELO): $(P_OBJ_VELO)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_POS): $(P_OBJ_POS)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_ATTITUDE): $(P_OBJ_ATTITUDE)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_PAYLOAD): $(P_OBJ_LOAD)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_GIF_3d): $(3D_GIF_OBJ)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_GIF_2d): $(2D_GIF_OBJ)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm -f $@ $^

clean:
	rm -f $(OBJS)
	rm -f DATA/*

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $< 
