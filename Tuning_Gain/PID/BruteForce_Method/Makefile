CC=gcc
CFLAGS=-Wall
HEADER=my_header.h
Z_DATA_DIR=Z_DATA/gain*
DATA_DIR=DATA/gain*
GARBAGE_F=DATA/gain*\ *
Z_GARBAGE_F=Z_DATA/gain*\ *

Z_OBJS=z_gain.o trajectries.o utils.o
CAL_OBJS=main.o trajectries.o utils.o
PLT_OBJ1=DRAW_FILE/plot_acc_v.o
PLT_OBJ2=DRAW_FILE/plot_acc_p.o
PLT_OBJ3=DRAW_FILE/plot_z.o

Z_GAIN=z_gain
CALC=calc
PLOT1=p_acc_v
PLOT2=p_acc_p
PLOT3=p_z

$(Z_GAIN): $(Z_OBJS)
	rm -f $(Z_DATA_DIR)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@
	rm -f $(Z_GARBAGE_F)

$(CALC):$(CAL_OBJS)
	rm -f $(DATA_DIR)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@
	rm -f $(GARBAGE_F)

$(PLOT1):$(PLT_OBJ1)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@ $^

$(PLOT2):$(PLT_OBJ2)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@ $^

$(PLOT3): $(PLT_OBJ3)
	$(CC) $(CFLAGS) -o $@ $^
	./$@
	rm $@ $^

%.o: %.c $(HEADER)
	$(CC) $(CFLAGS) -c -o $@ $<
