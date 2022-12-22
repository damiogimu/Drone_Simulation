# **************************************************************************** #
#                                                                              #
#                                                         :::      ::::::::    #
#    Makefile                                           :+:      :+:    :+:    #
#                                                     +:+ +:+         +:+      #
#    By: dmogi <mogi.t218140@gmail.com>             +#+  +:+       +#+         #
#                                                 +#+#+#+#+#+   +#+            #
#    Created: 2022/12/22 17:59:18 by dmogi             #+#    #+#              #
#    Updated: 2022/12/22 17:59:19 by dmogi            ###   ########.jp        #
#                                                                              #
# **************************************************************************** #

CC=gcc
CFLAGS=-Wall

C_FILES=sim_mba_pid_control.c utils.c
OBJS=sim_mba_pid_control.o utils.o

PLT_FILE1=DRAW_FILES/draw_error.c
PLT_FILE2=DRAW_FILES/draw_acc.c
PLT_FILE3=DRAW_FILES/draw_data.c
P_OBJ1=DRAW_FILES/draw_error.o
P_OBJ2=DRAW_FILES/draw_acc.o
P_OBJ3=DRAW_FILES/draw_data.o

CALC=calc
PLOT1=lerror # position error
PLOT2=lacc
PLOT_ATTITUDE=latti
PLOT_PAYLOAD=lpayload

$(CALC): $(C_FILES)
	$(CC) $(CFLAGS) -o $(CALC) $(C_FILES)
	./$(CALC)
	rm -f $(CALC)

$(PLOT1): $(P_OBJ1)
	rm -f $(P_OBJ1)
	./$(PLOT1)
	rm -f $(PLOT1)

$(PLOT2): $(P_OBJ2)
	rm -f $(P_OBJ2)
	./$(PLOT2)
	rm -f $(PLOT2)

$(PLOT_ATTITUDE): $(P_OBJ3)
	$(CC) $(CFLAGS) -D ATTI_F=1 -o $(PLOT_ATTITUDE) $(PLT_FILE3)
	rm -f $(P_OBJ3)
	./$(PLOT_ATTITUDE)
	rm -f $(PLOT_ATTITUDE)

$(PLOT_PAYLOAD): $(P_OBJ3)
	rm -f $(P_OBJ3)
	./$(PLOT_PAYLOAD)
	rm -f $(PLOT_PAYLOAD)

%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@
