MAKE   = make 
TARGET = my_code 
SOURCE = game_class.cpp
default: 
	mpicxx -o $(TARGET) $(SOURCE) --std=c++11
