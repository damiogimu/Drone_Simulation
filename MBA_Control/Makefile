CC=gcc
CFLAGS=-Wall

PLOT_DIR=Disp_files
DATA_DIR=Data
GRAPH_DIR=Graphs

HEADER=sim.h
OBJS=main.o trajectories.o utils.o
PLOT_OBJ=$(PLOT_DIR)/plot_data.o $(PLOT_DIR)/get_next_line.o $(PLOT_DIR)/get_next_line_utils.o
3D_GIF_OBJ=$(PLOT_DIR)/3d_anime.o
2D_GIF_OBJ=$(PLOT_DIR)/2d_anime.o

calc: $(OBJS)
	$(CC) -o $@ $^ -lm
	./$@
	rm $@

plot: $(PLOT_OBJ)
	$(CC) -o $@ $^
	./$@
	rm $@ $^

3d: $(3D_GIF_OBJ)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

2d: $(2D_GIF_OBJ)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

clean:
	rm -f $(OBJS)

fclean: clean
	rm -f $(DATA_DIR)/*
	rm -f $(GRAPH_DIR)/*

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $< 
