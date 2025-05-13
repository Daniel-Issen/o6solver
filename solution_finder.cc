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

#include "solution_finder.h"
#include "basis_consistency.h"
#include "pairing.h"
#include "constants.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <assert.h>

// Function to validate that a solution satisfies all clauses in the formula
bool validate_solution(const SATSolution& solution, 
                       const std::vector<std::vector<Literal>>& cnf_clauses) {
  std::cout << "Validating solution against the original problem..." << std::endl;
  
  // Check each clause in the formula
  for (size_t i = 0; i < cnf_clauses.size(); i++) {
    const auto& clause = cnf_clauses[i];
    bool clause_satisfied = false;
    
    // Check if at least one literal in the clause is satisfied
    for (const auto& literal : clause) {
      size_t var_idx = literal.var - 1;  // Convert to 0-indexed
      
      // Check if variable assignment satisfies this literal
      bool satisfied = false;
      if (var_idx < solution.assignments.size()) {
        // A positive literal (x_i) is satisfied if the variable is assigned true (1)
        // A negative literal (¬x_i) is satisfied if the variable is assigned false (-1)
        if ((solution.assignments[var_idx] == 1 && !literal.negated) ||
            (solution.assignments[var_idx] == -1 && literal.negated)) {
          satisfied = true;
        }
      }
      
      if (satisfied) {
        clause_satisfied = true;
        break;  // Only one literal needs to be true to satisfy the clause
      }
    }
    
    // If no literal satisfies this clause, the solution is invalid
    if (!clause_satisfied) {
      std::cout << "Validation failed: Clause " << (i + 1) << " is not satisfied." << std::endl;
      
      // Print the failing clause for debugging
      std::cout << "Failing clause: (";
      for (size_t j = 0; j < clause.size(); j++) {
        std::cout << clause[j].to_string();
        if (j < clause.size() - 1) std::cout << " ∨ ";
      }
      std::cout << ")" << std::endl;
      
      return false;
    }
  }
  
  std::cout << "Validation successful: All clauses are satisfied." << std::endl;
  return true;
}

SATSolution determine_solution
(std::vector<uint8_t>& basis_states,
 std::vector<uint8_t>& pair_states,
 std::vector<uint8_t>& term_states,
 Index n,
 int num_workers) {

  std::cout << "Attempting to determine a solution..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
    
  Index i = 0;
  Index j = i + 1;
  Index k = j + 1;
  
  // set three terms at a time
  Index starting_position = 0;
  for(i = 0; i < term_states.size() - 3; i += 3) {
    j = i + 1;
    k = j + 1;
    Index basis_idx = pair3d(i,j,k);
    uint8_t current_state = basis_states[basis_idx];
    // if the basis is either already fixed or trivially fixed,
    // continue to the next.
    update_basis_states(i, j, k,
			basis_idx,
			term_states, 
			pair_states, 
			basis_states);
    if(0 == (current_state & (current_state -1))) {
      continue;			// only one bit is set we can skip to
				// the next basis.
    }
    // pick the first valid solution
    basis_states[basis_idx] =
      basis_states[basis_idx] & -basis_states[basis_idx];

    // do a quick pass of updating each basis
    bool changed = false;
    do {
      changed = false;
      for(Index basis_index = starting_position;
	  basis_index < basis_states.size();
	  ++basis_index) {
	std::tuple<Index,Index,Index> ijk = unpair3d(basis_index);
	UpdateResult result =
	  update_basis_states(std::get<0>(ijk),
			      std::get<1>(ijk),
			      std::get<2>(ijk),
			      basis_index,
			      term_states,
			      pair_states,
			      basis_states);
	if(result.changed) {
	  changed = true;
	}
      }
    } while(changed);

    Index starting_basis_pair =
      pair2d(starting_position,starting_position + 1);
    Index ending_basis_pair =
      calculate_array_size_2d(calculate_array_size_3d(n));
    bool has_contradiction = false;
    if(num_workers < 2) {
      ensure_global_consistency(term_states,
				pair_states,
				basis_states,
				has_contradiction,
				starting_basis_pair,
				ending_basis_pair);
    } else {
      parallel_ensure_global_consistency(term_states,
					 pair_states,
					 basis_states,
					 has_contradiction,
					 starting_basis_pair,
					 ending_basis_pair,
					 num_workers);
    }
    starting_position = basis_idx;
  }
  SATSolution solution;
  solution.assignments.resize(n, 0);  // Initialize all as unassigned
  // we may have 1-2 terms still unset.
  j = i + 1;
  if(j < term_states.size()) { // two terms unset
    // update the pair based on its terms
    update_pair_states(i,j,term_states,pair_states);
    Index pair_idx = pair2d(i,j);
    uint8_t current_state = pair_states[pair_idx];
    if(!current_state) {
      std::cout << "this shouldn't happen" << std::endl;
      exit(0);
    }
    // pick the first valid solution
    pair_states[pair_idx] =
      pair_states[pair_idx] & -pair_states[pair_idx];
    update_pair_states(i,j,term_states,pair_states);
  } else if(i < term_states.size()) {
    if(term_states[i] == SET_ANY) { // one term unset
      term_states[i] = SET_POS;
    } else if(0 == term_states[i]) {
      std::cout << "this shouldn't happen" << std::endl;
    }
  }
  // Extract the solution from term_states
  for (Index i = 0; i < n; i++) {
    uint8_t state = term_states[i];
        
    if (state == SET_NEG) {
      solution.assignments[i] = -1;  // Negative assignment
    } else if (state == SET_POS) {
      solution.assignments[i] = 1;   // Positive assignment
    } else if (state == SET_ANY) {
      std::cout << "this shouldn't happen" << std::endl;
      exit(0);
    } else {
      std::cout << "this shouldn't happen" << std::endl;
      exit(0);
    }
  }
  
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Solution determination completed in "
	    << duration.count() << " ms" << std::endl;
    
  return solution;
}

void print_solution(const SATSolution& solution) {
  std::cout << "Solution:" << std::endl;
  for (size_t i = 0; i < solution.assignments.size(); i++) {
    // 1-indexed for output to match CNF file format
    switch (solution.assignments[i]) {
    case -1:
      std::cout << "0";
      break;
    case 1:
      std::cout << "1";
      break;
    default:
      std::cout << "UNDEFINED";
      break;
    }
    if(i < solution.assignments.size() -1) {
      std::cout << ",";
    }
  }
  std::cout << std::endl;
}

bool save_solution_to_file(const SATSolution& solution, 
			   const std::string& filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return false;
  }
    
  file << "# SAT problem solution" << std::endl;
  file << "# Variable assignments (1-indexed)" << std::endl;
    
  for (size_t i = 0; i < solution.assignments.size(); i++) {
    // Output in DIMACS-like format: positive or negative integers
    if (solution.assignments[i] == 1) {
      file << (i + 1);
    } else if (solution.assignments[i] == -1) {
      file << "-" << (i + 1);
    } else {
      file << "# x" << (i + 1) << " is undefined";
    }
    file << std::endl;
  }
    
  file << "0" << std::endl;  // End of solution marker
  return true;
}
