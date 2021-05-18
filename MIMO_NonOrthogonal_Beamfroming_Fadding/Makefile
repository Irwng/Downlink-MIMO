cc = g++
prom = main
deps = $(shell find ./ -name "*.h")
src = $(shell find ./ -name "*.cpp")
obj = $(src:%.cpp=%.o) 
 
$(prom): $(obj)
	$(cc) -static -m64 -o $(prom) $(obj)
%.o: %.cpp $(deps)
	$(cc) -Wall -static -O3 -DNDEBUG -m64 -c $< -I /usr/include/eigen3 -o $@
clean:
	rm -rf $(obj) $(prom)