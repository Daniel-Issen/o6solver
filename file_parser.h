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

#pragma once

#include <string>
#include <vector>
#include <random>

// Structure to represent a CNF literal
struct Literal {
    int var;     // Variable number (1-indexed as in DIMACS format)
    bool negated; // True if the literal is negated

    Literal(int v, bool neg) : var(v), negated(neg) {}
    
    std::string to_string() const {
        return (negated ? "-" : "") + std::to_string(var);
    }
};

// Parse a DIMACS CNF file
std::vector<std::vector<Literal>> parse_cnf_file(const std::string& filename, int& num_vars, int& num_clauses);

// Generate a random CNF formula
std::vector<std::vector<Literal>> generate_random_cnf(int num_vars, int num_clauses, int max_literals_per_clause, double negation_prob = 0.5);

// Overload that takes an explicit random generator
std::vector<std::vector<Literal>> generate_random_cnf(int num_vars, int num_clauses, int max_literals_per_clause, double negation_prob, std::mt19937& gen);
