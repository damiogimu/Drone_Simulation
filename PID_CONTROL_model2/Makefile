CC=gcc

FUNC_FILE1=sim_pid_control.c

PLT_FILE1=DRAW_FILES/draw_line_path.c
PLT_FILE2=DRAW_FILES/draw_error.c
PLT_FILE3=DRAW_FILES/draw_acc.c
PLT_FILE4=DRAW_FILES/draw_data.c

F_OBJ1=sim_pid_control.o

P_OBJ1=DRAW_FILES/draw_line_path.o
P_OBJ2=DRAW_FILES/draw_error.o
P_OBJ3=DRAW_FILES/draw_acc.o
P_OBJ4=DRAW_FILES/draw_data.o

LINE=line
ANIME=anime
PLOT=plot

PLOT1=lpath
PLOT2=lerror
PLOT3=lacc

DATA1=DATA/result_line_path
DATA2=DATA/result_line_error

ANIME_DIR=ANIMATION_PID
DATA_ANIME_PATH=../$(ANIME_DIR)/DATA/path_DATA
DATA_ANIME_CABLE=../$(ANIME_DIR)/DATA/cable_DATA
DATA_ANIME_XROTOR=../$(ANIME_DIR)/DATA/xrotor_DATA
DATA_ANIME_YROTOR=../$(ANIME_DIR)/DATA/yrotor_DATA

DATA_test_line=DATA/test_line

$(LINE): $(F_OBJ1)
	$(CC) $(CFLAGS) -o $(LINE) $(FUNC_FILE1) -lm
	./$(LINE) > $(DATA1)
	rm -f $(LINE)

test: $(F_OBJ1)
	$(CC) $(CFLAGS) -D CASE=1 -D TEST=1 -o $(LINE) $(FUNC_FILE1) -lm
	./$(LINE) > $(DATA_test_line)
	rm -f $(LINE)

$(PLOT1): $(P_OBJ1)
	$(CC) $(CFLAGS) -o $(PLOT1) $(PLT_FILE1)
	rm -f $(P_OBJ1)
	./$(PLOT1)
	rm -f $(PLOT1)

error: $(F_OBJ1)
	$(CC) $(CFLAGS) -D ERROR=1 -o $(LINE) $(FUNC_FILE1) -lm
	./$(LINE) > $(DATA2)
	rm -f $(LINE)

$(PLOT2): $(P_OBJ2)
	$(CC) $(CFLAGS) -o $(PLOT2) $(PLT_FILE2)
	rm -f $(P_OBJ2)
	./$(PLOT2)
	rm -f $(PLOT2)

$(PLOT3): $(P_OBJ3)
	$(CC) $(CFLAGS) -o $(PLOT3) $(PLT_FILE3)
	rm -f $(P_OBJ3)
	./$(PLOT3)
	rm -f $(PLOT3)

$(PLOT): $(P_OBJ4)
	$(CC) $(CFLAGS) -o $(PLOT) $(PLT_FILE4)
	rm -f $(P_OBJ4)
	./$(PLOT)
	rm -f $(PLOT)


$(ANIME): $(F_OBJ1)
	$(CC) $(CFLAGS) -D ANIME_F=1 -o $(ANIME) $(FUNC_FILE1) -lm
	./$(ANIME)
	rm -f $(ANIME)
