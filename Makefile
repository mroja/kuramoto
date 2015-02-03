
# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

CC = g++

TARGET = kuramoto_simulation
OBJS = kuramoto_simulation.o

$(TARGET): $(OBJS)
	$(CC) $^ -L$(BOOST_LIB) -o $@

$(OBJS): %.o: %.cpp
	$(CC) -Wall -Wextra -O3 -I$(BOOST_INC) -std=gnu++11 -c $< -o $@

clean:
	rm -f $(TARGET)
	rm -f *.o
	rm -f *.so

.PHONY: clean

