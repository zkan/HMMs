#
#	Makefile for compiling all code in this directory.
#	Author: Kan Ouivirach
#	Date: 17/08/2010
# 
CFLAGS = -g -Wall

#INC = -I /usr/include/opencv
#LIB = -lcxcore -lcv -lhighgui -lcvaux -lml
#LDFLAG = $(LIB)
#CC = g++ $(CFLAGS) $(INC)

CC = g++ $(CFLAGS)

#TARGET = ./hmm
#SRCS = hmm.cc baumwelch.cc forward.cc backward.cc viterbi.cc testCont.cc matrix.h

#PROG = genseq runForward runViterbi initHMM trainModel testCont runBaumWelch runBaumWelchWithSS onlineLearning runEnTraining runIncBaumWelch mergeHMMs inc-experiment testReadHMM test_get_files

PROG = test_get_files genseq run_forward run_viterbi init_hmm test_read_hmm run_baumwelch train_bw train_iml

all :	$(PROG)

# program to test the function that gets the list of files under a specific folder
test_get_files: test_get_files.o dirent.o 
	$(CC) test_get_files.o dirent.o -o test_get_files 

test_get_files.o: test_get_files.cc utility.h
	$(CC) -c test_get_files.cc

dirent.o: dirent.cc utility.h
	$(CC) -c dirent.cc
	
hmm.o: hmm.cc hmm.h
	$(CC) -c hmm.cc 

baumwelch.o: baumwelch.cc dirent.cc hmm.cc utility.h
	$(CC) -c baumwelch.cc

# program to generate a new sequence given a model
genseq: genseq.o hmm.o
	$(CC) genseq.o hmm.o -o genseq

# program to run the forward algorithm given a model
run_forward: run_forward.o hmm.o forward.o
	$(CC) run_forward.o hmm.o forward.o -o run_forward 

# program to run the Viterbi algorithm
run_viterbi: run_viterbi.o hmm.o viterbi.o
	$(CC) run_viterbi.o hmm.o viterbi.o -o run_viterbi

# program to initialize a new HMM
init_hmm: init_hmm.o hmm.o
	$(CC) init_hmm.o hmm.o -o init_hmm 

# program to test the function that reads a HMM file
test_read_hmm: test_read_hmm.o hmm.o
	$(CC) test_read_hmm.o hmm.o -o test_read_hmm 
	
# program to run the Baum-Welch algorithm to train only one sequence
run_baumwelch: run_baumwelch.o hmm.o baumwelch.o forward.o backward.o dirent.o
	$(CC) run_baumwelch.o hmm.o baumwelch.o forward.o backward.o dirent.o -o run_baumwelch

# program to run the Baum-Welch (BW) algorithm to train the whole data set
train_bw: train_bw.o hmm.o baumwelch.o forward.o backward.o dirent.o
	$(CC) train_bw.o hmm.o baumwelch.o forward.o backward.o dirent.o -o train_bw 
	
# program to run the incremental maximum-likelihood (IML) algorithm
train_iml: train_iml.o hmm.o iml.o forward.o backward.o dirent.o utility.h
	$(CC) train_iml.o hmm.o iml.o forward.o backward.o dirent.o -o train_iml

#testCont: testCont.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o testCont testCont.o hmm.o baumwelch.o forward.o backward.o -lm

#onlineLearning: onlineLearning.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o onlineLearning onlineLearning.o hmm.o baumwelch.o forward.o backward.o -lm

# The Ensemble Training algorithm
#runEnTraining: runEnTraining.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o runEnTraining runEnTraining.o hmm.o baumwelch.o forward.o backward.o -lm

# The Incremental Baum-Welch algorithm
#runIncBaumWelch: runIncBaumWelch.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o runIncBaumWelch runIncBaumWelch.o hmm.o baumwelch.o forward.o backward.o -lm

# Merge HMMs
#mergeHMMs: mergeHMMs.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o mergeHMMs mergeHMMs.o hmm.o baumwelch.o forward.o backward.o -lm

#inc-experiment: inc-experiment.o hmm.o baumwelch.o forward.o backward.o
#	$(CC) -o inc-experiment inc-experiment.o hmm.o baumwelch.o forward.o backward.o -lm


#test:
#	$(TARGET) model/normal.hmm sequence/normal.50.seq

clean:
	rm *.o $(PROG)







