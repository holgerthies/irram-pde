prefix=/Users/holgerthies/libraries/iRRAM/installed
exec_prefix=/Users/holgerthies/libraries/iRRAM/installed

CPPFLAGS             = -std=c++14 -O2 -Wall
LIBS = -liRRAM -liRRAMx -lmpfr -lgmp -lpng -L/Users/holgerthies/libraries/iRRAM/installed/lib -L/Users/holgerthies/work/iRRAMx/lib
INCLUDES = -I./src -I/Users/holgerthies/work/iRRAMx/include -I/Users/holgerthies/libraries/iRRAM/installed/include

CC := clang++ # This is the main compiler
BUILDDIR := build
TESTDIR := test
 
SRCEXT := cc
CFLAGS := -g -Wall -std=c++14 -Xlinker -rpath -Xlinker /Users/holgerthies/iRRAM/installed/lib


$(TESTDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(TESTDIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -lboost_unit_test_framework -o $@ $< 

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

