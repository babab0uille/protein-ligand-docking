#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cassert>

using namespace std;

const int GRID_SIZE = 10;         // Grid size
const int POPULATION_SIZE = 100;  // Number of ligands in each generation
const int GENERATIONS = 50;       // Number of generations to run
const double MUTATION_RATE = 0.1; // Probability of mutation

// Structure to represent a ligand
struct Ligand
{
    int x, y;       // Position of the ligand on the grid
    double fitness; // Fitness score
};

// Function to calculate fitness of a ligand
double calculateFitness(const Ligand &ligand, const vector<pair<int, int>> &activeSite)
{
    double score = 0.0;
    for (const auto &site : activeSite)
    {
        double distance = sqrt(pow(ligand.x - site.first, 2) + pow(ligand.y - site.second, 2));
        score += 1.0 / (1.0 + distance); // Higher score for closer distances
    }
    return score;
}

// Initialize population with random ligand positions
vector<Ligand> initializePopulation(int gridSize, int populationSize)
{
    vector<Ligand> population;
    for (int i = 0; i < populationSize; i++)
    {
        Ligand ligand = {rand() % gridSize, rand() % gridSize, 0.0};
        population.push_back(ligand);
    }
    return population;
}

// Selection: Choose the top individuals based on fitness
vector<Ligand> selectFittest(const vector<Ligand> &population, int numToSelect)
{
    vector<Ligand> selected(population.begin(), population.begin() + numToSelect);
    return selected;
}

// Crossover: Combine two parent ligands to create a child ligand
Ligand crossover(const Ligand &parent1, const Ligand &parent2)
{
    Ligand child;
    child.x = (parent1.x + parent2.x) / 2;
    child.y = (parent1.y + parent2.y) / 2;
    return child;
}

// Mutation: Introduce small random changes
void mutate(Ligand &ligand, int gridSize)
{
    if ((rand() / (double)RAND_MAX) < MUTATION_RATE)
    {
        ligand.x += rand() % 3 - 1; // Random shift [-1, 1]
        ligand.y += rand() % 3 - 1; // Random shift [-1, 1]

        // Bound check
        ligand.x = max(0, min(gridSize - 1, ligand.x));
        ligand.y = max(0, min(gridSize - 1, ligand.y));
    }
}

// Function to export generation data to CSV
void exportResultsToCSV(const vector<Ligand> &population, int generation)
{
    ofstream file("results.csv", ios::app);
    if (generation == 1)
    { // Write the header for the first generation
        file << "Generation,X,Y,Fitness\n";
    }
    for (const auto &ligand : population)
    {
        file << generation << "," << ligand.x << "," << ligand.y << "," << ligand.fitness << "\n";
    }
    file.close();
}

// Main Genetic Algorithm function
Ligand geneticAlgorithm(const vector<pair<int, int>> &activeSite)
{
    vector<Ligand> population = initializePopulation(GRID_SIZE, POPULATION_SIZE);

    for (int gen = 0; gen < GENERATIONS; ++gen)
    {
        for (auto &ligand : population)
        {
            ligand.fitness = calculateFitness(ligand, activeSite);
        }
        sort(population.begin(), population.end(), [](const Ligand &a, const Ligand &b)
             { return a.fitness > b.fitness; });
        vector<Ligand> newPopulation = selectFittest(population, POPULATION_SIZE / 2);

        while (newPopulation.size() < POPULATION_SIZE)
        {
            const Ligand &parent1 = newPopulation[rand() % (POPULATION_SIZE / 2)];
            const Ligand &parent2 = newPopulation[rand() % (POPULATION_SIZE / 2)];
            Ligand child = crossover(parent1, parent2);
            mutate(child, GRID_SIZE);
            newPopulation.push_back(child);
        }

        population = newPopulation;

        cout << "Generation " << gen + 1 << " - Best Ligand: ("
             << population[0].x << ", " << population[0].y
             << ") with fitness " << population[0].fitness << endl;

        exportResultsToCSV(population, gen + 1);
    }
    return population[0];
}

// Test Cases
void runTests()
{
    // Test case 1
    vector<pair<int, int>> testActiveSite1 = {{5, 5}, {7, 7}, {3, 3}};
    Ligand result1 = geneticAlgorithm(testActiveSite1);
    assert(result1.x >= 0 && result1.x < GRID_SIZE);
    assert(result1.y >= 0 && result1.y < GRID_SIZE);

    // Test case 2
    vector<pair<int, int>> testActiveSite2 = {{1, 1}, {2, 2}, {3, 3}};
    Ligand result2 = geneticAlgorithm(testActiveSite2);
    assert(result2.x >= 0 && result2.x < GRID_SIZE);
    assert(result2.y >= 0 && result2.y < GRID_SIZE);

    // Test case 3
    vector<pair<int, int>> testActiveSite3 = {{2, 3}, {4, 5}, {6, 7}};
    Ligand result3 = geneticAlgorithm(testActiveSite3);
    assert(result3.x >= 0 && result3.x < GRID_SIZE);
    assert(result3.y >= 0 && result3.y < GRID_SIZE);

    // Test case 4
    vector<pair<int, int>> testActiveSite4 = {{8, 8}, {0, 0}, {3, 3}};
    Ligand result4 = geneticAlgorithm(testActiveSite4);
    assert(result4.x >= 0 && result4.x < GRID_SIZE);
    assert(result4.y >= 0 && result4.y < GRID_SIZE);

    // Test case 5
    vector<pair<int, int>> testActiveSite5 = {{4, 5}, {6, 6}, {9, 9}};
    Ligand result5 = geneticAlgorithm(testActiveSite5);
    assert(result5.x >= 0 && result5.x < GRID_SIZE);
    assert(result5.y >= 0 && result5.y < GRID_SIZE);

    cout << "All tests passed!" << endl;
}

// Main Function
int main()
{
    srand(static_cast<unsigned int>(time(0))); // Seed for randomness

    // Active site points
    vector<pair<int, int>> activeSite = {{5, 5}, {6, 6}, {4, 4}};

    cout << "Running Genetic Algorithm for Protein-Ligand Docking...\n";

    Ligand bestLigand = geneticAlgorithm(activeSite);

    cout << "\nOptimal Ligand Position: (" << bestLigand.x << ", " << bestLigand.y
         << ") with fitness " << bestLigand.fitness << endl;

    // Run test cases
    cout << "\nRunning Test Cases...\n";
    runTests();

    return 0;
}
