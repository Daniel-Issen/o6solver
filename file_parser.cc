// MIT License

// Copyright (c) 2025 Daniel Issen

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "file_parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <random>
#include <algorithm>

// Parse a CNF file (DIMACS or similar format)
std::vector<std::vector<Literal>> parse_cnf_file(const std::string& filename,
						 int& num_vars,
						 int& num_clauses) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + filename);
  }
    
  std::vector<std::vector<Literal>> clauses;
  std::string line;
  int max_var_id = 0;
    
  // Read through the file
  while (std::getline(file, line)) {
    // Skip comments and empty lines
    if (line.empty() || line[0] == 'c') {
      continue;
    }
        
    // Handle problem line if present (for backward compatibility)
    if (line[0] == 'p') {
      std::istringstream iss(line);
      std::string p, cnf;
      int declared_vars, declared_clauses;
      iss >> p >> cnf >> declared_vars >> declared_clauses;
            
      // Store these values but don't rely on them
      std::cout << "Found problem line: " << declared_vars << " variables, " 
		<< declared_clauses << " clauses declared" << std::endl;
      continue;
    }
        
    // Parse clause line
    std::istringstream iss(line);
    int var;
    std::vector<Literal> clause;
        
    while (iss >> var && var != 0) {
      // Track the maximum variable ID to determine the total number
      // of variables
      int var_id = std::abs(var);
      max_var_id = std::max(max_var_id, var_id);
            
      // DIMACS format: positive numbers represent positive literals,
      // negative numbers represent negated literals
      clause.push_back(Literal(var_id, var < 0));
    }
        
    if (!clause.empty()) {
      clauses.push_back(clause);
    }
  }
    
  // Set the derived values
  num_vars = max_var_id;
  num_clauses = static_cast<int>(clauses.size());
    
  std::cout << "Derived from input: " << num_vars << " variables, " 
	    << num_clauses << " clauses" << std::endl;
    
  return clauses;
}

// Generate a random CNF formula with explicit random generator
std::vector<std::vector<Literal>> generate_random_cnf
(int num_vars, 
 int num_clauses, 
 int max_literals_per_clause, 
 double negation_prob,
 std::mt19937& gen) {
  
  std::uniform_int_distribution<int> var_dist(1, num_vars);
  std::uniform_int_distribution<int> size_dist(1, max_literals_per_clause);
  std::uniform_real_distribution<double> neg_dist(0.0, 1.0);
    
  std::vector<std::vector<Literal>> clauses;
    
  for (int i = 0; i < num_clauses; i++) {
    int clause_size = size_dist(gen);
    std::vector<Literal> clause;
    std::set<int> used_vars; // To avoid duplicate variables in the same clause
        
    for (int j = 0; j < clause_size; j++) {
      int var;
      // Find an unused variable for this clause
      do {
	var = var_dist(gen);
      } while (used_vars.find(var) != used_vars.end());
            
      used_vars.insert(var);
      bool negated = (neg_dist(gen) < negation_prob);
      clause.push_back(Literal(var, negated));
    }
        
    clauses.push_back(clause);
  }
    
  return clauses;
}

// Maintain the original signature for backward compatibility
std::vector<std::vector<Literal>> generate_random_cnf
(int num_vars, 
 int num_clauses, 
 int max_literals_per_clause, 
 double negation_prob) {
  std::random_device rd;
  std::mt19937 gen(rd());
  return generate_random_cnf(num_vars, 
			     num_clauses, 
			     max_literals_per_clause,
			     negation_prob,
			     gen);
}
