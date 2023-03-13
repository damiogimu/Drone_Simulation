CC=gcc
CFLAGS=-Wall

PLOT_DIR=GNUPLOT_FILES

HEADER=my_header.h
OBJS=main.o utils.o
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
	rm -f DATA/*
	rm -f GIF/*

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
