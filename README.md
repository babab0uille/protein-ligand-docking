# Protein-Ligand Docking Genetic Algorithm

This repository implements a genetic algorithm to optimize ligand positions for protein-ligand docking. The algorithm simulates evolution by selecting, crossing over, and mutating candidate ligand positions based on their proximity to active sites, aiming to find the optimal docking position.

## Overview

The genetic algorithm follows the principles of natural selection:
1. **Initialization**: A random population of ligands (representing potential docking positions) is generated.
2. **Fitness Calculation**: Each ligand is evaluated based on its proximity to known active sites. A fitness score is computed, with closer ligands having a higher score.
3. **Selection**: The top-performing ligands are selected for reproduction.
4. **Crossover**: Pairs of ligands combine to form new ligands, mimicking genetic recombination.
5. **Mutation**: Random changes are introduced to simulate genetic variation.
6. **Iteration**: The process continues for a set number of generations to converge on the optimal solution.

## How to Run the Code

### Prerequisites
Ensure you have the following installed:
- C++ compiler (e.g., `g++`)
- CMake (optional, for building)

### To compile and run the algorithm:
1. Clone the repository:
    ```bash
    git clone https://github.com/babab0uille/protein-ligand-docking.git
    cd protein-ligand-docking
    ```

2. Compile the code:
    ```bash
    g++ -o docking_algorithm ProteinLigandDocking.cpp
    ```

3. Run the algorithm:
    ```bash
    ./docking_algorithm
    ```
This will execute the genetic algorithm on predefined active sites and output the results of each generation. It will also run internal test cases to validate the algorithm's correctness.

## How the Tests Work

The tests are embedded in the `main` function and will automatically run when you execute the program. These tests validate that the algorithm produces ligands within the defined grid size, ensuring correctness. If all tests pass, you will see: " All tests passed! "

## Saving Results to a CSV File

The results of each generation, including the positions of the ligands and their fitness scores, are saved to a CSV file named `results.csv`. This file can be opened in Excel or any spreadsheet software for easy analysis. 
