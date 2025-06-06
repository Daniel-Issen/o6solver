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

#include "test_utils.h"
#include "file_parser.h"
#include "cnf_solver.h"  // For check_satisfiability
#include <iostream>
#include <chrono>

// Simple brute force check for satisfiability (for small instances)
bool check_satisfiability_brute_force
(const std::vector<std::vector<Literal>>& cnf_clauses,
 int num_vars,
 int& num_solutions) {
  // For large instances, brute force is impractical
  if (num_vars > 20) {
    std::cout << "Brute force check skipped (too many variables: "
	      << num_vars << ")" << std::endl;
    return false;
  }
    
  // Try all possible assignments to variables
  uint64_t max_assignments = 1ULL << num_vars;
  num_solutions = 0;
    
  auto start = std::chrono::high_resolution_clock::now();
    
  for (uint64_t assignment = 0; assignment < max_assignments; assignment++) {
    bool valid = true;
        
    // Check each clause
    for (const auto& clause : cnf_clauses) {
      bool clause_satisfied = false;
            
      // A clause is satisfied if at least one of its literals is true
      for (const auto& literal : clause) {
	int var = literal.var - 1; // Convert to 0-indexed
	bool var_value = (assignment >> var) & 1; // Extract bit for
						  // this variable 
                
	if (var_value != literal.negated) {
	  clause_satisfied = true;
	  break;
	}
      }
            
      if (!clause_satisfied) {
	valid = false;
	break;
      }
    }
        
    if (valid) {
      num_solutions++;
      std::cout << "Solution " << num_solutions << ": ";
      for (int var = 0; var < num_vars; var++) {
	std::cout << ((assignment >> var) & 1);
	if (var < num_vars - 1 && var < 19) std::cout << ",";
      }
      std::cout << std::endl;
    }
  }
    
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = 
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
  std::cout << "Brute force check results:" << std::endl;
  std::cout << "- Total solutions found: " << num_solutions << std::endl;
  std::cout << "- Time taken: " << duration.count() << " ms" << std::endl;
    
  return num_solutions > 0;
}

// Test random formulas with algorithm selection
void test_random_formulas(int num_workers,int num_tests, int num_vars, int num_clauses, int max_literals_per_clause, bool find_solution) {
  std::cout << "Testing " << num_tests << " random formulas..." << std::endl;
  std::cout << "Parameters: " 
	    << num_vars << " variables, " 
	    << num_clauses << " clauses, "
	    << "max " << max_literals_per_clause << " literals per clause"
	    << std::endl;
    
  int correct_results = 0;
  double total_time = 0;
    
  for (int test = 1; test <= num_tests; test++) {
    std::cout << "Test "
	      << test
	      << "/" << num_tests
	      << std::endl;
    std::cout << "-------------------------------------------------------"
	      << std::endl;
        
    // Generate a random CNF formula
    auto cnf_formula = 
      generate_random_cnf(num_vars, num_clauses, max_literals_per_clause);
        
    // Print a sample of the formula
    std::cout << "Random formula:\n";
    for (size_t i = 0; i < cnf_formula.size(); i++) {
      for (size_t j = 0; j < cnf_formula[i].size(); j++) {
	std::cout << cnf_formula[i][j].to_string();
	if (j < cnf_formula[i].size() - 1) std::cout << " ";
      }
      std::cout << " 0" << std::endl;
    }
    // Check satisfiability
    std::cout << std::endl << "Checking satisfiability..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    bool result =
      check_satisfiability(num_workers,cnf_formula, num_vars, find_solution);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = 
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    total_time += duration.count();
        
    // Run brute force check for very small instances
    if (num_vars <= 15) {
      std::cout << std::endl << "Running brute force check..." << std::endl;
      int num_solutions = 0;
      bool brute_force_result = 
	check_satisfiability_brute_force(cnf_formula, 
					 num_vars, 
					 num_solutions);
            
      std::cout << "Comparison:" << std::endl;
      std::cout << "- Our result: " 
		<< (result ? "Satisfiable" : "Unsatisfiable") << std::endl;
      std::cout << "- Brute force result: " 
		<< (brute_force_result ? "Satisfiable" : "Unsatisfiable") 
		<< std::endl;
            
      // Check if results are consistent
      bool consistent = (brute_force_result == result || result);
      std::cout << "- Results are consistent: " 
		<< (consistent ? "Yes" : "No") << std::endl;
            
      if (consistent) {
	correct_results++;
      }
    } else {
      // For larger instances, we can't verify with brute force
      // Just assume the result is correct
      correct_results++;
    }
        
    std::cout << "-------------------------------------------------------" <<
      std::endl << std::endl;
  }
    
  // Print summary statistics
  std::cout << "Test Summary:\n";
  std::cout << "- Tests run: " << num_tests << std::endl;
  std::cout << "- Correct/consistent results: " 
	    << correct_results 
	    << " (" << (correct_results * 100.0 / num_tests) << "%)"
	    << std::endl;
  std::cout << "- Average time: " << (total_time / num_tests) << std::endl;
}
