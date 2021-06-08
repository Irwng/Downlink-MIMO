cc = g++
prom = main
deps = $(shell find ./ -name "*.h")
src = $(shell find ./ -name "*.cpp")
obj = $(src:%.cpp=%.o) 
 
$(prom): $(obj)
	$(cc) -static -m64 -o $(prom) $(obj)
%.o: %.cpp $(deps)
	$(cc) -Wall -O3 -DNDEBUG -m64 -c $< -I /usr/include/eigen3 -o $@
clean:
	rm -rf $(obj) $(prom)
	
# main: main.o header.o viterbi.o map.o
# 	g++ main.o header.o viterbi.o map.o -o main
# main.o: main.cpp header.h
# 	g++ -c main.cpp -I C:\\Users\\wangyi\\eigen3 -o main.o
# header.o: header.cpp header.h
# 	g++ -c header.cpp -I C:\\Users\\wangyi\\eigen3 -o header.o
# viterbi.o: viterbi.cpp header.h
# 	g++ -c viterbi.cpp -I C:\\Users\\wangyi\\eigen3 -o viterbi.o
# map.o: map.cpp header.h
# 	g++ -c map.cpp -I C:\\Users\\wangyi\\eigen3 -o map.o


# /usr/include/eigen3
# .PHONY:clean
# clean:
#     rm -rf *.o test
