%.o: %.cpp
	g++ -c $<

main: main.o UserClass.o Checkpoint.o
	g++ -o main main.o UserClass.o Checkpoint.o

clean:
	rm *.o
