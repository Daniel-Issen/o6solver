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

#include "parallel_solver.h"
#include <iostream>
#include <chrono>
#include <algorithm>
#include "cnf_solver.h"
#include "solution_finder.h"

// Function to divide work among workers
std::vector<WorkSegment> divide_work(uint64_t n, int num_workers) {
  std::vector<WorkSegment> segments;
    
  // Calculate total number of bases for logging
  uint64_t total_bases = calculate_array_size_3d(n);
  std::cout << "Total basis triplets: " << total_bases << std::endl;
    
  // Use a balanced division strategy - divide the total basis count evenly
  uint64_t bases_per_worker = (total_bases + num_workers - 1) / num_workers;
  std::cout << "Bases per worker (target): " << bases_per_worker << std::endl;
    
  uint64_t current_position = 0;
    
  for (int worker = 0; worker < num_workers; worker++) {
    WorkSegment segment;
    segment.start_position = current_position;
        
    // Calculate end position for this worker
    uint64_t end_position =
      std::min(current_position + bases_per_worker, total_bases);
        
    // Convert positions to (i,j,k) coordinates for debugging
    uint64_t start_i, start_j, start_k;
    uint64_t end_i, end_j, end_k;
        
    if (current_position < total_bases) {
      std::tie(start_i, start_j, start_k) = unpair3d(current_position);
    } else {
      // This worker gets no work
      start_i = start_j = start_k = 0;
    }
        
    if (end_position > 0 && end_position <= total_bases) {
      std::tie(end_i, end_j, end_k) = unpair3d(end_position - 1);
    } else {
      end_i = n - 2;
      end_j = n - 1;
      end_k = n;
    }
        
    // Make sure the last worker gets all remaining work
    if (worker == num_workers - 1) {
      end_i = n - 2;
      end_j = n - 1;
      end_k = n;
    }
        
    segment.max_i = end_i;
    segment.max_j = end_j;
    segment.max_k = end_k;
        
    // Update current position for next worker
    current_position = end_position;
        
    // Log worker assignments for debugging
    std::cout << "Worker " << worker << ": "
	      << "Start pos: " << segment.start_position << " "
	      << "(" << start_i << "," << start_j << "," << start_k << ") -> "
	      << "Max: (" << segment.max_i << "," << segment.max_j << ","
	      << segment.max_k << ") " << "Bases: "
	      << (end_position - segment.start_position) << std::endl;
        
    segments.push_back(segment);
  }
    
  return segments;
}

// Worker function that processes a segment
WorkerResult process_segment(const WorkSegment& segment,
			     uint64_t n,
			     const std::vector<uint8_t>& term_states,
			     const std::vector<uint8_t>& pair_states,
			     const std::vector<uint8_t>& basis_states) {
    
  WorkerResult result;
  result.term_states = term_states;
  result.pair_states = pair_states;
  result.basis_states = basis_states;
  result.has_contradiction = false;
    
  // Process this segment
  ensure_global_consistency(n,
			    result.term_states,
			    result.pair_states,
			    result.basis_states,
			    result.has_contradiction,
			    segment.start_position,
			    segment.max_i,
			    segment.max_j,
			    segment.max_k);
  return result;
}

// Function to merge worker results
bool merge_worker_results(
			  std::vector<WorkerResult>& worker_results,
			  std::vector<uint8_t>& term_states,
			  std::vector<uint8_t>& pair_states,
			  std::vector<uint8_t>& basis_states,
			  bool& has_contradiction) {
    
  bool changed = false;
  has_contradiction = false;  // Initialize to false
    
  // Initialize with the first worker's results
  term_states = worker_results[0].term_states;
  pair_states = worker_results[0].pair_states;
  basis_states = worker_results[0].basis_states;
    
  if (worker_results[0].has_contradiction) {
    has_contradiction = true;
    return false; // No need to continue if contradiction found
  }
    
  // Merge results from other workers
  for (size_t i = 1; i < worker_results.size(); i++) {
    auto& worker = worker_results[i];
        
    if (worker.has_contradiction) {
      has_contradiction = true;
      return false; // No need to continue if contradiction found
    }
        
    // Intersection of allowed states (bitwise AND)
    for (size_t j = 0; j < term_states.size(); j++) {
      uint8_t old_term = term_states[j];
      term_states[j] &= worker.term_states[j];
      if (term_states[j] != old_term) changed = true;
      if (term_states[j] == 0) {
	has_contradiction = true;
	return false; // Contradiction detected during merge
      }
    }
        
    for (size_t j = 0; j < pair_states.size(); j++) {
      uint8_t old_pair = pair_states[j];
      pair_states[j] &= worker.pair_states[j];
      if (pair_states[j] != old_pair) changed = true;
      if (pair_states[j] == 0) {
	has_contradiction = true;
	return false; // Contradiction detected during merge
      }
    }
        
    for (size_t j = 0; j < basis_states.size(); j++) {
      uint8_t old_basis = basis_states[j];
      basis_states[j] &= worker.basis_states[j];
      if (basis_states[j] != old_basis) changed = true;
      if (basis_states[j] == 0) {
	has_contradiction = true;
	return false; // Contradiction detected during merge
      }
    }
  }
    
  return changed;
}

// Main parallel solver function
bool parallel_check_satisfiability
(const std::vector<std::vector<Literal>>& cnf_clauses,
 int num_vars,
 int num_workers,
 bool find_solution,
 const std::string& solution_file) {
    
  if (num_workers <= 1) {
    // Fallback to sequential solver if only one worker
    return check_satisfiability(cnf_clauses, 
				num_vars, 
				find_solution, 
				solution_file);
  }
    
  std::cout << "Using parallel solver with "
	    << num_workers << " workers" << std::endl;
    
  // Initialize state arrays
  std::vector<uint8_t> term_states(num_vars, SET_ANY);
  std::vector<uint8_t> pair_states(calculate_array_size_2d(num_vars),
				   SET_ANY_ANY);
  std::vector<uint8_t> basis_states(calculate_array_size_3d(num_vars),
				    SET_ANY_ANY_ANY);
    
  // Apply constraints directly
  int working_num_vars = num_vars;
  bool initial_consistency = apply_constraints(
					       cnf_clauses, 
					       working_num_vars, 
					       term_states, 
					       pair_states, 
					       basis_states
					       );
    
  if (!initial_consistency) {
    std::cout << "Formula is unsatisfiable (detected during initial constraint application)" << std::endl;
    return false;
  }
    
  // Parallel global consistency check
  auto start = std::chrono::high_resolution_clock::now();
    
  bool has_contradiction = false;
  bool changed = true;
  int iterations = 0;
    
  while (changed && !has_contradiction) {
    iterations++;
    std::cout << "Iteration " << iterations << "..." << std::endl;
        
    // Divide work among workers
    auto work_segments = divide_work(working_num_vars, num_workers);
        
    // Spawn worker threads
    std::vector<std::future<WorkerResult>> futures;
    for (const auto& segment : work_segments) {
      futures.push_back(std::async(std::launch::async, 
				   process_segment, 
				   segment, 
				   working_num_vars, 
				   term_states, 
				   pair_states, 
				   basis_states));
    }
        
    // Collect results
    std::vector<WorkerResult> worker_results;
    for (auto& future : futures) {
      worker_results.push_back(future.get());
    }
        
    // Merge results
    changed = merge_worker_results(worker_results, 
				   term_states, 
				   pair_states, 
				   basis_states, 
				   has_contradiction);
        
        
    if (worker_results.size() > 0 && worker_results[0].has_contradiction) {
      has_contradiction = true;
    }
        
  }
    
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
    std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
  std::cout << "Results:\n";
  std::cout << "- Iterations: " << iterations << std::endl;
  std::cout << "- Contradiction detected: "
	    << (has_contradiction ? "Yes" : "No") << std::endl;
  std::cout << "- Time taken: " << duration.count() << " ms" << std::endl;
    
  // If a contradiction was detected, the formula is unsatisfiable
  if (has_contradiction) {
    return false;
  }
    
  // If we want to find a solution and no contradiction was detected
  if (find_solution) {
    SATSolution solution = determine_solution(
					      basis_states, 
					      pair_states,
					      term_states,
					      num_vars
					      );
        
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
