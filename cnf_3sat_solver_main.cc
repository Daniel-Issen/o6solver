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

#include <iostream>
#include <string>
#include <stdexcept>
#include "file_parser.h"
#include "cnf_solver.h"
#include "test_utils.h"

#if 1

#include "parallel_solver.h"

// Modified main function in cnf_3sat_solver_main.cc
int main(int argc, char* argv[]) {
  std::cout << "Optimized CNF Solver\n";
  std::cout << "===================\n\n";
    
  // Parse command line arguments
  bool run_tests = false;
  bool find_solution = false;
  std::string cnf_file;
  std::string solution_file;
  int num_tests = 10;
  int test_vars = 10;
  int test_clauses = 20;
  int max_literals = 3;
  int num_workers = 1;  // Default to sequential execution
    
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--test" || arg == "-t") {
      run_tests = true;
      if (i + 1 < argc && argv[i+1][0] != '-') {
        num_tests = std::stoi(argv[i+1]);
        i++;
      }
    } else if (arg == "--vars" || arg == "-v") {
      if (i + 1 < argc) {
        test_vars = std::stoi(argv[i+1]);
        i++;
      }
    } else if (arg == "--clauses" || arg == "-c") {
      if (i + 1 < argc) {
        test_clauses = std::stoi(argv[i+1]);
        i++;
      }
    } else if (arg == "--literals" || arg == "-l") {
      if (i + 1 < argc) {
        max_literals = std::stoi(argv[i+1]);
        i++;
      }
    } else if (arg == "--solve" || arg == "-s") {
      find_solution = true;
    } else if (arg == "--output" || arg == "-o") {
      if (i + 1 < argc) {
        solution_file = argv[i+1];
        i++;
      }
    } else if (arg == "--workers" || arg == "-w") {  // New parameter
      if (i + 1 < argc) {
        num_workers = std::stoi(argv[i+1]);
        i++;
      }
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " [options] [cnf_file]\n";
      std::cout << "Options:\n";
      std::cout << "  --test, -t [num]       Run tests on random formulas (default: 10 tests)\n";
      std::cout << "  --vars, -v [num]       Number of variables for random tests (default: 10)\n";
      std::cout << "  --clauses, -c [num]    Number of clauses for random tests (default: 20)\n";
      std::cout << "  --literals, -l [num]   Maximum literals per clause (default: 3)\n";
      std::cout << "  --solve, -s            Find and output a solution if formula is satisfiable\n";
      std::cout << "  --output, -o [file]    Save solution to the specified file\n";
      std::cout << "  --workers, -w [num]    Number of worker threads for parallel execution (default: 1)\n";
      std::cout << "  --help, -h             Show this help message\n";
      return 0;
    } else {
      cnf_file = arg;
    }
  }
    
  try {
    if (run_tests) {
      // Run tests on random formulas
      test_random_formulas(num_tests,
                           test_vars,
                           test_clauses,
                           max_literals,
                           find_solution);
    } else if (!cnf_file.empty()) {
      // Parse the CNF file
      int num_vars, num_clauses;
      std::cout << "Parsing CNF file: " << cnf_file << std::endl;
      auto cnf_clauses = parse_cnf_file(cnf_file, num_vars, num_clauses);
            
      std::cout << "Formula details:\n";
      std::cout << "- Variables: " << num_vars << std::endl;
      std::cout << "- Clauses: " << num_clauses << std::endl;
            
      // Print a sample of the clauses
      std::cout << "\nSample clauses:\n";
      for (size_t i = 0; i < std::min(cnf_clauses.size(), size_t(10)); i++) {
        std::cout << "(";
        for (size_t j = 0; j < cnf_clauses[i].size(); j++) {
          std::cout << cnf_clauses[i][j].to_string();
          if (j < cnf_clauses[i].size() - 1) std::cout << " ∨ ";
        }
        std::cout << ")\n";
      }
      if (cnf_clauses.size() > 10) {
        std::cout << "... and " << (cnf_clauses.size() - 10)
                  << " more clauses\n";
      }
            
      // Check satisfiability
      std::cout << "\nChecking satisfiability...\n";
      bool result =
        check_satisfiability(num_workers,
			     cnf_clauses,
			     num_vars,
			     find_solution,
			     solution_file);
            
      // Run brute force check for small instances
      int num_solutions = 0;
      bool brute_force_result = false;
      if (num_vars <= 20) {
        std::cout << "\nRunning brute force check...\n";
        brute_force_result = check_satisfiability_brute_force(
							      cnf_clauses,
							      num_vars,
							      num_solutions
							      );
      }
            
      // Report final result
      std::cout << "\nFinal result:\n";
      if (num_vars <= 20 && !result && !brute_force_result) {
        std::cout << "Formula is UNSATISFIABLE (confirmed by brute force)\n";
      } else if (num_vars <= 20 && brute_force_result) {
        std::cout << "Formula is SATISFIABLE with "
                  << num_solutions << " solutions\n";
      } else if (!result) {
        std::cout << "Formula is UNSATISFIABLE\n";
      } else {
        std::cout << "Formula is SATISFIABLE\n";
                
        // If we're asked to find a solution but haven't already done so
        if (find_solution && !result) {
          std::cout << "Cannot determine a solution as the formula is unsatisfiable.\n";
        }
      }
    } else {
      // If no file specified and no test requested, run a small demo
      std::cout << "No CNF file specified. Running a small demo with random formula...\n\n";
      test_random_formulas(3, 8, 15, 3, find_solution);
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
    
  return 0;
}

#else

int main(int argc, char* argv[]) {
  std::cout << "Optimized CNF Solver\n";
  std::cout << "===================\n\n";
    
  // Parse command line arguments
  bool run_tests = false;
  bool find_solution = false;
  std::string cnf_file;
  std::string solution_file;
  std::string augmented_file;
  int num_tests = 10;
  int test_vars = 10;
  int test_clauses = 20;
  int max_literals = 3;
    
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];
    if (arg == "--test" || arg == "-t") {
      run_tests = true;
      if (i + 1 < argc && argv[i+1][0] != '-') {
	num_tests = std::stoi(argv[i+1]);
	i++;
      }
    } else if (arg == "--vars" || arg == "-v") {
      if (i + 1 < argc) {
	test_vars = std::stoi(argv[i+1]);
	i++;
      }
    } else if (arg == "--clauses" || arg == "-c") {
      if (i + 1 < argc) {
	test_clauses = std::stoi(argv[i+1]);
	i++;
      }
    } else if (arg == "--literals" || arg == "-l") {
      if (i + 1 < argc) {
	max_literals = std::stoi(argv[i+1]);
	i++;
      }
    } else if (arg == "--solve" || arg == "-s") {
      find_solution = true;
    } else if (arg == "--output" || arg == "-o") {
      if (i + 1 < argc) {
	solution_file = argv[i+1];
	i++;
      }
    } else if (arg == "--augmented" || arg == "-a") {
      if (i + 1 < argc) {
	augmented_file = argv[i+1];
	i++;
      }
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " [options] [cnf_file]\n";
      std::cout << "Options:\n";
      std::cout << "  --test, -t [num]       Run tests on random formulas (default: 10 tests)\n";
      std::cout << "  --vars, -v [num]       Number of variables for random tests (default: 10)\n";
      std::cout << "  --clauses, -c [num]    Number of clauses for random tests (default: 20)\n";
      std::cout << "  --literals, -l [num]   Maximum literals per clause (default: 3)\n";
      std::cout << "  --solve, -s            Find and output a solution if formula is satisfiable\n";
      std::cout << "  --output, -o [file]    Save solution to the specified file\n";
      std::cout << "  --augmented, -a [file] Generate augmented CNF file with additional constraints\n";
      std::cout << "  --help, -h             Show this help message\n";
      return 0;
    } else {
      cnf_file = arg;
    }
  }
    
  try {
    if (run_tests) {
      // Run tests on random formulas
      test_random_formulas(num_tests,
			   test_vars,
			   test_clauses,
			   max_literals,
			   find_solution);
    } else if (!cnf_file.empty()) {
      // Parse the CNF file
      int num_vars, num_clauses;
      std::cout << "Parsing CNF file: " << cnf_file << std::endl;
      auto cnf_clauses = parse_cnf_file(cnf_file, num_vars, num_clauses);
            
      std::cout << "Formula details:\n";
      std::cout << "- Variables: " << num_vars << std::endl;
      std::cout << "- Clauses: " << num_clauses << std::endl;
            
      // Print a sample of the clauses
      std::cout << "\nSample clauses:\n";
      for (size_t i = 0; i < std::min(cnf_clauses.size(), size_t(10)); i++) {
	std::cout << "(";
	for (size_t j = 0; j < cnf_clauses[i].size(); j++) {
	  std::cout << cnf_clauses[i][j].to_string();
	  if (j < cnf_clauses[i].size() - 1) std::cout << " ∨ ";
	}
	std::cout << ")\n";
      }
      if (cnf_clauses.size() > 10) {
	std::cout << "... and " << (cnf_clauses.size() - 10)
		  << " more clauses\n";
      }
            
      // Check satisfiability
      std::cout << "\nChecking satisfiability...\n";
      bool result =
	check_satisfiability(cnf_clauses,
			     num_vars,
			     find_solution,
			     solution_file, 
			     !augmented_file.empty() ? cnf_file : "",
			     augmented_file);
            
      // Run brute force check for small instances
      int num_solutions = 0;
      bool brute_force_result = false;
      if (num_vars <= 20) {
	std::cout << "\nRunning brute force check...\n";
	brute_force_result =
	  check_satisfiability_brute_force(cnf_clauses,
					   num_vars,
					   num_solutions);
      }
            
      // Report final result
      std::cout << "\nFinal result:\n";
      if (num_vars <= 20 && !result && !brute_force_result) {
	std::cout << "Formula is UNSATISFIABLE (confirmed by brute force)\n";
      } else if (num_vars <= 20 && brute_force_result) {
	std::cout << "Formula is SATISFIABLE with "
		  << num_solutions << " solutions\n";
      } else if (!result) {
	std::cout << "Formula is UNSATISFIABLE\n";
      } else {
	std::cout << "Formula is SATISFIABLE\n";
                
	// If we're asked to find a solution but haven't already done so
	if (find_solution && !result) {
	  std::cout <<
	    "Cannot determine a solution as the formula is unsatisfiable.\n";
	}
      }
            
      // If an augmented file was generated, print info about it
      if (!augmented_file.empty() && result) {
	std::cout << "\nAugmented CNF file has been generated: "
		  << augmented_file << std::endl;
	std::cout << "This file contains additional clauses derived from the consistency algorithm." << std::endl;
	std::cout << "You can use this augmented file with any SAT solver for potentially faster solving." << std::endl;
      }
    } else {
      // If no file specified and no test requested, run a small demo
      std::cout << "No CNF file specified. Running a small demo with random formula...\n\n";
      test_random_formulas(3, 8, 15, 3, find_solution);
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
    
  return 0;
}
#endif
