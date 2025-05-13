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

#include "cnf_solver.h"
#include "constants.h"
#include "pairing.h"
#include "basis_consistency.h"
#include "solution_finder.h"
#include "parallel_solver.h"
#include <chrono>
#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>

// Directly apply CNF constraints without creating unnecessary dummy variables
bool apply_constraints(const std::vector<std::vector<Literal>>& cnf_clauses,
		       int& num_vars,
		       std::vector<uint8_t>& term_states,
		       std::vector<uint8_t>& pair_states,
		       std::vector<uint8_t>& basis_states) {  
  uint64_t max_var_id = num_vars;
  for(const auto& clause : cnf_clauses) {
    std::vector<Literal> sorted_clause = clause;
    auto comp = [](const Literal &a, const Literal &b) {
      return (a.var < b.var);
    };
    std::sort(sorted_clause.begin(),sorted_clause.end(),comp);
    if(sorted_clause.size() == 0) continue;
    if(sorted_clause.size() == 1) {
      term_states[sorted_clause[0].var -1] &=
	oned_clear_masks[sorted_clause[0].negated];
      if(!(term_states[sorted_clause[0].var -1])) {
	return false;
      }
    } else if(sorted_clause.size() == 2) {
      uint64_t idx = pair2d(sorted_clause[0].var -1,
			    sorted_clause[1].var -1);
      pair_states[idx] &=
	twod_clear_masks[sorted_clause[0].negated][sorted_clause[1].negated];
      if(!(pair_states[idx])) {
	return false;
      }
    } else if(sorted_clause.size() == 3) {
      uint64_t idx = pair3d(sorted_clause[0].var -1,
			    sorted_clause[1].var -1,
			    sorted_clause[2].var -1);
      basis_states[idx] &=
	threed_clear_masks
	[sorted_clause[0].negated]
	[sorted_clause[1].negated]
	[sorted_clause[2].negated];
      if(!(basis_states[idx])) {
	return false;
      }
    } else {
      // Handle clauses with more than 3 literals
      // We'll break them down into multiple 3-literal clauses with
      // auxiliary variables 
      // For a clause (a ∨ b ∨ c ∨ d ∨ e ∨ ...)
      // We introduce auxiliary variables z1, z2, ... and create:
      // (a ∨ b ∨ z1), (¬z1 ∨ c ∨ z2),
      // (¬z2 ∨ d ∨ z3), ..., (¬z_n ∨ last ∨ second_last) 

      // resize arrays
      term_states.resize((term_states.size() + sorted_clause.size() - 3),
			 SET_ANY);
      pair_states.
	resize(calculate_array_size_2d(term_states.size()),SET_ANY_ANY);
      basis_states.
	resize(calculate_array_size_3d(term_states.size()),SET_ANY_ANY_ANY);
      // Handle first clause (a ∨ b ∨ z1)
      uint64_t idx = pair3d(sorted_clause[0].var -1,
			    sorted_clause[1].var -1,
			    max_var_id);
      
      basis_states[idx] &=
	threed_clear_masks
	[(int)sorted_clause[0].negated]
	[(int)sorted_clause[1].negated]
	[0];
      
      max_var_id++; // Increment after use
      
      // Handle intermediate clauses (¬z_i ∨ term ∨ z_{i+1})
      for (size_t i = 2; i < sorted_clause.size() - 2; i++) {
        uint64_t prev_var_id = max_var_id - 1;
        idx = pair3d(sorted_clause[i].var -1,
		     prev_var_id,
		     max_var_id);
        basis_states[idx] &=
          threed_clear_masks
          [(int)sorted_clause[i].negated]
          [1]
          [0];
        max_var_id++;
      }
      
      // Handle last clause ( -z_i v (last_term -1) v (last_term)
      uint64_t prev_var_id = max_var_id - 1;
      idx = pair3d(sorted_clause[sorted_clause.size() -2].var -1,
		   sorted_clause[sorted_clause.size() -1].var -1,
		   prev_var_id);
      basis_states[idx] &=
	threed_clear_masks
	[sorted_clause[sorted_clause.size() -2].negated]
	[sorted_clause[sorted_clause.size() -1].negated]
	[1];
    }
  }
  num_vars = max_var_id;
  return true; // No contradictions found during initial constraint application
}

// Check satisfiability
bool check_satisfiability
(int num_workers,
 const std::vector<std::vector<Literal>>& cnf_clauses, 
 int num_vars, 
 bool find_solution,
 const std::string& solution_file) {
  // Initialize state arrays
  std::vector<uint8_t> term_states(num_vars, SET_ANY);
  std::vector<uint8_t> pair_states(calculate_array_size_2d(num_vars), 
				   SET_ANY_ANY);
  std::vector<uint8_t> basis_states(calculate_array_size_3d(num_vars), 
				    SET_ANY_ANY_ANY);
  // Apply constraints directly
  int working_num_vars = num_vars;

  bool initial_consistency =
    apply_constraints(cnf_clauses, 
		      working_num_vars, 
		      term_states, 
		      pair_states, 
		      basis_states);
    
  if (!initial_consistency) {
    std::cout << "Formula is unsatisfiable (detected during initial constraint application)" << std::endl;
    return false;
  }
  // Run global consistency check
  bool has_contradiction = false;
  uint64_t ending_basis_pair =
    calculate_array_size_2d(calculate_array_size_3d(working_num_vars));
  auto start = std::chrono::high_resolution_clock::now();
  if(num_workers < 2) {
    ensure_global_consistency(term_states, 
			      pair_states, 
			      basis_states, 
			      has_contradiction,
			      0,ending_basis_pair);
  } else {
    parallel_ensure_global_consistency(term_states, 
				       pair_states, 
				       basis_states, 
				       has_contradiction,
				       0,ending_basis_pair,
				       num_workers);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = 
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
  std::cout << "Results:\n";
  std::cout << "- Contradiction detected: "
	    << (has_contradiction ? "Yes" : "No") << std::endl;
  std::cout << "- Time taken: " << duration.count() << " ms" << std::endl;
    
  // If a contradiction was detected, the formula is unsatisfiable
  if (has_contradiction) {
    return false;
  }
  // If we want to find a solution and no contradiction was detected
  if (find_solution) {
    SATSolution solution = 
      determine_solution(basis_states, 
			 pair_states,
			 term_states,
			 num_vars,
			 0);
			 //			 num_workers);

    // Validate the solution against the original problem
    bool valid = validate_solution(solution, cnf_clauses);
    if (valid) {
      std::cout << "verified solution" << std::endl;
    } else {
      std::cerr << "Warning: The determined solution does not satisfy the formula!" << std::endl;
      // This shouldn't happen if the algorithm is correct
      return false;
    }
    
    // Print the solution to console
    print_solution(solution);
        
    // Save to file if requested
    if (!solution_file.empty()) {
      if (save_solution_to_file(solution, solution_file)) {
	std::cout << "Solution saved to file: " << solution_file << std::endl;
      } else {
	std::cerr << "Failed to save solution to file." << std::endl;
      }
    }
  }
    
  // If no contradiction was found, the formula is satisfiable
  return true;
}
