CC = gcc
CUNIT = /usr/local/Cellar/cunit/2.1-2
CUNITLIB = -I$(CUNIT)/include -L$(CUNIT)/lib 
OBJ = no
MYMATH = mymath
MYLINESEARCH = line_search

$(OBJ):
	$(CC) $(CUNITLIB) test_$(OBJ).c ../src/$(OBJ).c -lcunit -o test_$(OBJ).o

$(MYMATH):
	$(CC) $(CUNITLIB) test_$(MYMATH).c ../src/$(MYMATH).c -lcunit -o test_$(MYMATH).o

$(MYLINESEARCH):
	$(CC) $(CUNITLIB) test_$(MYLINESEARCH).c ../src/$(MYLINESEARCH).c ../src/$(MYMATH).c -lcunit -o test_$(MYLINESEARCH).o

clean:
	rm -f *.o

