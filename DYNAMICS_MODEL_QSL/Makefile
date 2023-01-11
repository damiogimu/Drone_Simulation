CC=gcc
CFLAG=-Wall
GNL=TOOL/get_next_line.c TOOL/get_next_line_utils.c 
RESULT_LAG=./RESULT/result_lag
RESULT_NEWTON=./RESULT/result_newton
RESULT_ERROR=./RESULT/result_error

LAG_OBJ=sim_lagrange.o
NEWTON_OBJ=sim_newton.o
ERROR_OBJ=error.o
PLOT_OBJ=draw.o
PLOT_FILE=draw.c

OBJS=$(LAG_OBJ) $(NEWTON_OBJ) $(ERROR_OBJ) $(PLOT_OBJ)

CALC_LAG=calc_lag
CALC_NEWTON=calc_newton
CALC_ERROR=calc_error

PLOT_LAG_POS=plot_lag
PLOT_NEWTON_POS=plot_newton
PLOT_ERROR_POS=plot_error
PLOT_ACC=plot_acc


$(CALC_NEWTON):$(NEWTON_OBJ)
	$(CC) $(CFLAG) -o $(CALC_NEWTON) $(NEWTON_OBJ)
	./$(CALC_NEWTON) > $(RESULT_NEWTON)
	rm $(CALC_NEWTON)	

$(CALC_LAG):$(LAG_OBJ)
	$(CC) $(CFLAG) -o $(CALC_LAG) $(LAG_OBJ)
	./$(CALC_LAG) > $(RESULT_LAG)
	rm $(CALC_LAG)

$(CALC_ERROR):$(ERROR_OBJ)
	$(CC) $(CFLAG) -o $(CALC_ERROR) $(ERROR_OBJ) $(GNL)
	./$(CALC_ERROR) > $(RESULT_ERROR)
	rm $(CALC_ERROR)

$(PLOT_LAG_POS):$(PLOT_OBJ)
	$(CC) $(CFLAG) -D LAG_F=1 -o $(PLOT_LAG_POS) $(PLOT_FILE)
	./$(PLOT_LAG_POS)
	rm $(PLOT_LAG_POS)

$(PLOT_NEWTON_POS):$(PLOT_OBJ)
	$(CC) $(CFLAG) -o $(PLOT_NEWTON_POS) $(PLOT_FILE)
	./$(PLOT_NEWTON_POS)
	rm $(PLOT_NEWTON_POS)

$(PLOT_ERROR_POS):$(PLOT_OBJ)
	$(CC) $(CFLAG) -D ERROR_F=1 -o $(PLOT_ERROR_POS) $(PLOT_FILE)
	./$(PLOT_ERROR_POS)
	rm $(PLOT_ERROR_POS)

$(PLOT_ACC):$(PLOT_OBJ)
	$(CC) $(CFLAG) -D ACC_F=1 -o $(PLOT_ACC) $(PLOT_FILE)
	./$(PLOT_ACC)
	rm $(PLOT_ACC)

clean:$(OBJS)
	rm $^

%.o:%.c
	$(CC) $(CFLAG) -c $<
