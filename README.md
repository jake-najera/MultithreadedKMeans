# Multi-threaded K-Means
Taking Naive K-Means and applying three methods of multithreading, with the intention to compare runtime between each method for CS 355.
1. Single-threaded
2. Multi-threaded Centroid Update (MTCU)
3. Multi-threaded Workers Parallelization (MTWP) 

# Single-threaded K-Means
1. Assign samples to a centroid
2. Update the centroids coordinates to center of samples assigned.
# Multi-threaded Centroid Update
This version multithreads the second step of K-Means.
1. Assign samples to a centroid.
2. Update centroid coordinates on threads.
3. A barrier, to keep threads in sync. 


![image](https://user-images.githubusercontent.com/54923430/112529221-09329f00-8d73-11eb-9867-3705c295d7e9.png)


Note: This method suffers from parent thread bottlenecking causing runtime performance to not be optimal.

# Multi-threaded Workers Parallelization
This version implements paralellism by divvying both steps in algorithm to workers.
1. Delegates slices of the samples array to T worker threads, to assign to centroids.
2. A barrier, to keep threads in sync.
3. Delegates slices of centroids array to T worker threads, to update coordinates.
4. A barrier, to keep threads in sync.


![image](https://user-images.githubusercontent.com/54923430/112529278-1780bb00-8d73-11eb-8eed-6200b9b2d990.png)

# Usage
Make with provided makefile.

Five/Six Arguments:
1. An integer specifying the type of threading:
	+ 0 - Single
	+ 1 - Parallel Centroid Update (MTCU)
	+ 2 - Concurrent-Means (MTWP)

+ If and only if type == 2, supply number of threads, T,  to use

2. Number of Clusters, K
3. Number of Samples, N
4. Minimum of Random Values When Generating Samples On Fixed Seed
5. Maximum of Random Values When Generating Samples On Fixed Seed

Example Usage using a Bash prompt:

Single threaded, 2 Clusters, 500 Samples, [-1000, 1000]
```bash                
./multikmeans.exe 0 2 500 -1000 1000
```
MTCU, 4 Clusters, 500 Samples, [-1000, 1000]
```bash  
./multikmeans.exe 1 4 500 -1000 1000
```
MTWP, 6 Threads, 6 Clusters, 500 Samples, [-1000, 1000]
```bash  
./multikmeans.exe 2 6 6 500 -1000 1000
```
Example Data Collection Command Using Bash:
```bash
for j in {1000..20000..1000};
	do for i in {1..10};
		do time ./multikmeans.exe 2 10 10 ${j} -1000 1000; 
		echo ${j};
    done;
done
```
