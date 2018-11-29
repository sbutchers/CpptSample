
# If you already have your own Makefile you can 
# replace all of this, but you need to keep this line
# at the top:
include ${COSMOSIS_SRC_DIR}/config/compilers.mk


#For python, leave this as is.  For C/C++/Fortran, 
#comment the first line and uncomment the second:
#all: 
all: cppt_sample.so

# ... and then provide the command to compile it, perhaps
# like this for a simple example in C:
cppt_sample.so: cosmosis_module.cpp
	$(CXX) -std=c++14 -shared -I/home/sean/cosmosis/ -o cppt_sample.so cosmosis_module.cpp


#Replace this to put in some kind of test command for
#your module, if you have one.
test:
	@echo "Alas, cppt_module/cppt_sample module has no tests"

#Add anything else to clean here
clean:
	rm -f *.o *.mod

