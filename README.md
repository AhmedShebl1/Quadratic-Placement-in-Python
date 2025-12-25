# Quadratic-Placement-in-Python

Quadratic-Placement-in-Python is a Python implementation of a Recursive Quadratic Placement (QP) engine for VLSI physical design.  
It minimizes total wirelength by formulating placement as a system of linear equations and solving it using an analytical, force-directed approach.

The engine uses top-down recursive partitioning with pseudo-pad propagation to handle boundary constraints and preserve global connectivity during local optimization.

---

## Overview

This project implements an Analytical Quadratic Placer commonly used in VLSI physical design flows.  
It distributes gates across a chip area by solving quadratic wirelength minimization problems and recursively refining placement regions.

Key characteristics:
- Analytical (not simulated annealing)
- Force-directed quadratic formulation
- Recursive top-down placement
- Clique-based net modeling

---

## Core Features

### Quadratic Programming Formulation
- Minimizes squared wirelength
- Solves linear systems of the form:
  A x = b
- X and Y coordinates are solved independently

### Recursive Partitioning
- Top-down placement strategy
- Alternates between vertical and horizontal cuts
- Divides placement regions into equal-capacity subregions

### Clique Model Connectivity
- Multi-terminal nets are converted into weighted point-to-point edges
- Each edge weight is 1 / (k - 1), where k is the net degree
- Provides accurate quadratic wirelength estimation

---

## Algorithm Flow

The placement engine follows a standard analytical placement pipeline:

1. Global Placement  
   - Solve the quadratic system for all gates

2. Slicing  
   - Sort gates by position
   - Split the region into two equal-capacity subregions

3. Terminal Propagation  
   - External connections are converted into pseudo-pads
   - Pseudo-pads act as fixed anchors during local solves

4. Recursive Refinement  
   - Repeat until the specified recursion depth is reached

---

## Input and Output Formats

### Input File Format (.txt)

The input describes the circuit netlist and pad locations.

&lt;NumGates&gt; &lt;NumNets&gt;

&lt;GateID&gt; &lt;NumNetsConnected&gt; &lt;Net1&gt; &lt;Net2&gt;

...

&lt;NumPads&gt;

&lt;PadID&gt; &lt;NetConnected&gt; &lt;PadX&gt; &lt;PadY&gt;

...

### Output File Format

The output lists optimized coordinates for each gate.

&lt;GateID&gt; &lt;X_Coordinate&gt; &lt;Y_Coordinate&gt;

...

---

## Execution Guide

### Prerequisites

- Python 3.x
- NumPy
- SciPy (for sparse linear solvers)

### Installation

pip install numpy scipy

### Running the Placer

Run the engine by providing the input netlist and recursion depth:

python qp_engine.py <input_file_path> <depth>

Example:

python qp_engine.py benchmarks/ckt1.txt 4
---

## Project Scope

This project demonstrates key concepts in VLSI physical design, including:
- Analytical quadratic placement
- Wirelength minimization
- Sparse linear system solving
- Recursive top-down partitioning
- Pseudo-pad based constraint handling
