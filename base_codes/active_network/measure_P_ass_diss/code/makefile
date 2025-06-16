SRCS=	random.c main.c energy.c init.c criticalcasimir.c analysis.c terminate_block.c switchfunction.c quaternion.c storage.c version_specific_functions.c

OBJS=	random.o main.o energy.o init.o criticalcasimir.o analysis.o terminate_block.o switchfunction.o quaternion.o  storage.o version_specific_functions.o


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







