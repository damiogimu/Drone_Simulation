CC=gcc
CFLAGS=-Wall

HEADER=my_header.h
OBJS=main.o trajectories.o utils.o
PLOT_OBJ=GNUPLOT_FILES/plot_data.o GNUPLOT_FILES/get_next_line.o GNUPLOT_FILES/get_next_line_utils.o
3D_GIF_OBJ=GNUPLOT_FILES/3d_anime.o
2D_GIF_OBJ=GNUPLOT_FILES/2d_anime.o

CALC=calc
PLOT_DATA=pdf
PLOT_GIF_3d=3d
PLOT_GIF_2d=2d

$(CALC): $(OBJS)
	$(CC) -o $@ $^ -lm
	./$@
	rm $@

$(PLOT_DATA): $(PLOT_OBJ)
	$(CC) -o $@ $^
	./$@
	rm $@ $^	

$(PLOT_GIF_3d): $(3D_GIF_OBJ)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

$(PLOT_GIF_2d): $(2D_GIF_OBJ)
	$(CC) -o $@ $^
	./$@
	rm -f $@ $^

clean:
	rm -f $(OBJS)

%.o:%.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
