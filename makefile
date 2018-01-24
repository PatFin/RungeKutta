all : FirstOrder.exe SecondOrder.exe RobotArm.exe

%.exe : %.o gnuplot_i.o
	gcc -o $@ $^ -lm

%.o : %.c
	gcc -c -o $@ $<

clean :
	@rm *.exe
	@rm *.o
