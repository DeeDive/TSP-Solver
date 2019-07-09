# TSP-Solver

**Traveling Salesman Problem** is one of well-known **NP-Complete** problems in the field of computer science. 

In this project, we attempt to solve the problem by using **various classic algorithms** and **Qt 5.12 for visualization**.

---
**Environment**
+ C++ 11/14/17
+ Qt 5.12

---

### Overview 

<p align="center">
	<img src="https://github.com/Infinitestudy/TSP-Solver/blob/master/demo-img/overview.jpg" alt="Overview"  width="70%" height="70%">
	<p align="center">
		<em>Figure 1. Overview of TSP-Solver</em>
	</p>
</p>

---
### History ###

In the past few decades, many researchers have endeavored to solve the NPC problems, especially TSP. 

There has been a shifted focus from certain complexity improvement to specific problem instances.

As a consequence, the scale of the problem we can solve is growing exponentially (see Figure 2).


<p align="center">
	<img src="https://github.com/Infinitestudy/TSP-Solver/blob/master/demo-img/his.jpg" alt="History"  width="70%" height="70%">
	<p align="center">
		<em>Figure 2. History of TSP</em>
	</p>
</p>

---

### Data
TSP solver supports both **random generated data** and **standardized data from TSPLIB OR VLSI Data Sets** .

---
### Algorithms

We implement nine classic algorithms to solve the TSP problem. Four of them are algorithms for exact solution, and the other five are algorithms for approximate solution. From my perspective, there is a trade-off between precision and time for algorithms that give approximate solution.

---
### Demo for instance with more cities

<p align="center">
	<img src="https://github.com/Infinitestudy/TSP-Solver/blob/master/demo-img/morecities.jpg" alt="morecities_instance"  width="70%" height="70%">
	<p align="center">
		<em>Figure 3. Solve instance with more cities</em>
	</p>
</p>
