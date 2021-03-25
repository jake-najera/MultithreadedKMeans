# Macros ===============================================================

CC=g++
CFLAGS=-lpthread -std=c++11
EXE=multikmeans
FOLDER = 
ERASE=rm
FILES= $(FOLDER)Common.cpp $(FOLDER)driver.cpp $(FOLDER)K_meansClustering.cpp $(FOLDER)Centroid.cpp $(FOLDER)SampleUnsupervised.cpp $(FOLDER)MTCU_Clustering.cpp $(FOLDER)ParallelClustering.cpp
# Build target =========================================================

$(OUTDIR)$(EXE) :
	$(CC) -o $(EXE) $(FILES) $(CFLAGS)

# Other targets ========================================================
	
clean :
	$(ERASE) $(EXE)

rebuild :
	$(MAKE) clean
	$(MAKE)

