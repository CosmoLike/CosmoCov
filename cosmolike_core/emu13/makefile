P = main.o emu.o hubble.o emu_noh.o
          
LDFLAGS = -O

emu.exe: $(P) $(CSTUFF)
	gcc  -o emu.exe  $(LDFLAGS) $(P) -lgsl -lgslcblas -lm 

clean:
	/bin/rm -f *.o

