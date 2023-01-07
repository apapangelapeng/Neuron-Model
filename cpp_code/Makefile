# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -Wall -g
 
# ****************************************************
#exectuables
.PHONY : all
all :static_ap nullcline propogate_ap reduced_ap static_wt_ap static_vcompare





static_ap: static_ap.o
	$(CC) $(CFLAGS) -o static_ap static_ap.o
	"dont forget to pass a current value when calling the executable like ./static_ap 0.1"

nullcline:nullcline.o
	$(CC) $(CFLAGS) -o nullcline nullcline.o



propogate_ap: propogate_ap.o
	$(CC) $(CFLAGS) -o propogate_ap propogate_ap.o


reduced_ap: reduced_ap.o
	$(CC) $(CFLAGS) -o reduced_ap reduced_ap.o

static_wt_ap:static_wt_ap.o
	$(CC) $(CFLAGS) -o static_wt_ap static_wt_ap.o


static_ap_hcn:static_ap_hcn.o
	$(CC) $(CFLAGS) -o static_ap_hcn static_ap_hcn.o


static_vcompare:static_vcompare.o
	$(CC) $(CFLAGS) -o static_vcompare static_vcompare.o


static_hcn_ap:static_hcn_ap.o
	$(CC) $(CFLAGS) -o static_hcn_ap static_hcn_ap.o
# .o aka object files 
 


static_ap.o: Static_AP.cpp 
	$(CC) $(CFLAGS) Static_AP.cpp -c -o static_ap.o

propogate_ap.o: Propogate_AP.cpp 
	$(CC) $(CFLAGS) propogate_AP.cpp -c -o propogate_ap.o

reduced_ap.o: Reduced_AP.cpp
	$(CC) $(CFLAGS) Reduced_AP.cpp -c -o reduced_ap.o


nullcline.o: nullcline.cpp
	$(CC) $(CFLAGS) nullcline.cpp -c -o nullcline.o


static_wt_ap.o: static_wt_ap.cpp
	$(CC) $(CFLAGS) static_wt_ap.cpp -c -o static_wt_ap.o

static_ap_hcn.o: static_ap_hcn.cpp
	$(CC) $(CFLAGS) static_ap_hcn.cpp -c -o static_ap_hcn.o

static_vcompare.o: static_vcompare.cpp
	$(CC) $(CFLAGS) static_vcompare.cpp -c -o static_vcompare.o

static_hcn_ap.o: static_hcn_ap.cpp
	$(CC) $(CFLAGS) static_hcn_ap.cpp -c -o static_hcn_ap.o

.PHONY : clean
clean :
	-rm *.o $(objects) 
	rm -rf *.ll *.exe *.s *.o testall.log
	
.PHONY : clean_img
clean_img:
	rm -rf ../graphs/*.png

clean_csv:
	rm -rf -rf ../data_files/*.csv