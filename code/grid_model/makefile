CC          	= gcc
RM          	= rm
CFLAGS      	= -O3 -lm -fPIC
LDFLAGS         = -lm
SOURCES_MAIN    = Voltage_Control_main.c
SOURCES     	= Voltage_Control_interface.c sim_code/rtGetInf.c sim_code/rtGetNaN.c sim_code/rt_nonfinite.c sim_code/Voltage_Control.c sim_code/Voltage_Control_data.c 
OBJECT_MAIN		= Voltage_Control_main.o
OBJECTS     	= $(SOURCES:.c=.o)
EXECS       	= $(SOURCES:%.c=%)
TARGET			= Voltage_Control_main
SHARED_LIB		= Voltage_Control_lib.so
decimation		= 1
stop_time 		= -1


.PHONY: all
all: $(OBJECT_MAIN) $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECT_MAIN) $(OBJECTS) -o $(TARGET)

.PHONY: shared_lib
shared_lib: $(OBJECTS)
	$(CC) -fPIC -shared $(CFLAGS) $(OBJECTS) -o $(SHARED_LIB)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

.o.: 
	$(CC) $^ $(LDFLAGS) -o $@

run: all
	./$(TARGET) $(stop_time) $(decimation)
	
.PHONY: clean
clean:
	-@ $(RM) $(OBJECT_MAIN)
	-@ $(RM) $(OBJECTS) 
	-@ $(RM) $(TARGET)
	-@ $(RM) $(SHARED_LIB)