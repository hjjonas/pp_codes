SRCS=	 main.c init.c random.c criticalcasimir.c analysis.c quaternion.c switchfunction.c mc.c energy.c storage.c terminate_block.c bmd.c constructnn.c version_specific_functions.c
OBJS=	 main.o init.o random.o criticalcasimir.o analysis.o quaternion.o switchfunction.o mc.o energy.o storage.o terminate_block.h bmd.o constructnn.o version_specific_functions.o

HFILES = 

#

CF =   gcc
CF2 =   gcc-10
FFLAGS = -Wno-deprecated-declarations   
PARALLEL = -fopenmp
LIBS =   -lm    -framework OpenGL -framework Glut -framework Cocoa 
OUT =	../runfile.run


opt:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -o $(OUT) $(SRCS:.c=.o) $(LIBS)

prof:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -pg -o $(OUT) $(SRCS) $(LIBS)

opt3:	$(SRCS)
	$(CF) -Wno-deprecated-declarations  -O3 -o $(OUT) $(SRCS) $(LIBS) 

opt4:	$(SRCS)
	$(CF2) -Wno-deprecated-declarations  -O3 -o $(OUT) $(SRCS) $(LIBS) $(PARALLEL)

deb:	$(SRCS)
	$(CF) -g -O2 -o rods.deb $(SRCS) $(LIBS)

$(SRCS):	$(HFILES)
$(OBJS):	$(HFILES)	







