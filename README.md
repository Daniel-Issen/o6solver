# o6solver
an O(n^6) 3sat solver.

Usage: ./cnf_3sat_solver [options] [cnf_file]
Options:
  --test, -t [num]       Run tests on random formulas (default: 10 tests)
  --vars, -v [num]       Number of variables for random tests (default: 10)
  --clauses, -c [num]    Number of clauses for random tests (default: 20)
  --literals, -l [num]   Maximum literals per clause (default: 3)
  --solve, -s            Find and output a solution if formula is satisfiable
  --output, -o [file]    Save solution to the specified file
  --workers, -w [num]    Number of worker threads for parallel execution (default: 1)
  --help, -h             Show this help message

This program is hillariously slow but does run in polynomial time.  It
is consistently outperformed by minisat.  You won't be using this to
factor RSA keys anytime soon.  The algorithm is easily parallelized
however and this program implements a simple parallelization strategy.

To build, use the included Makefile.  The program relies on
__uint128_t, a non-standard data type implemented by g++ and clang.

enjoy

Daniel

Daniel Issen
daniel.issen@gmail.com
