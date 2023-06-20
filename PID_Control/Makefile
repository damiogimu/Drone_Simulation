CC=gcc
CFLAGS=-Wall

PFILE_DIR=Disp_files
DATA_DIR=Data
GRAPH_DIR=Graphs

HEADER=sim.h
OBJS=main.o trajectories.o utils.o
PLOT_OBJ=$(PFILE_DIR)/plot_data.o $(PFILE_DIR)/get_next_line.o $(PFILE_DIR)/get_next_line_utils.o
3D_GIF_OBJ=$(PFILE_DIR)/3d_anime.o
2D_GIF_OBJ=$(PFILE_DIR)/2d_anime.o

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
