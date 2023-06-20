CC=gcc
CFLAGS=-Wall
HEADER=my_header.h
DATA_DIR=DATA/gain*
GARBAGE_F=DATA/gain*\ *

CAL_OBJS=main.o trajectries.o utils.o
PLT_OBJ1=plot_acc_v.o
PLT_OBJ2=plot_acc_p.o

CALC=calc
PLOT1=p_acc_v
PLOT2=p_acc_p

$(CALC):$(CAL_OBJS)
	rm -f $(DATA_DIR)
	$(CC) $(CFALGS) -o $@ $^
	./$@
	rm $@
	rm -f $(GARBAGE_F)

$(PLOT1):$(PLT_OBJ1)
	$(CC) $(CFALGS) -o $@ $^
	./$@
	rm $@ $^

$(PLOT2):$(PLT_OBJ2)
	$(CC) $(CFALGS) -o $@ $^
	./$@
	rm $@ $^

%.o: %.c $(HEADER)
	$(CC) $(CFALGS) -c $<
