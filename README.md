# Matrix Game Designer

## Solving an instance

Given a matrix instance A, cost vectors c_r and c_s, and a budget B, the "solve_game" function returns the optimal strategy x and decision vectors r and s. If given vectors r_fix and s_fix, rows and columns corresponding to 1 (0) will be purchased (not purchased) in the strategy, respectively.

Hueristic approaches are implemented in the "solve_game_greedy_MILP" and "solve_game_greedy_frequency" functions. When conditionally dominated rows are known, they can be incorporated to strengthen the formulation via the "solve_game_with_cuts" function.

## Generating random instances

The "src/parameters.jl" file contains all functions for generating random matrix game instances. There are four "create_matrix" methods, each of which receives an entry range from which to draw entries from as well as a seed for reproducibility. You may pass either the number of rows and the number of columns, or create a square matrix by passing one such quantity. Both of these have an additional method to pass a numerical vector R so that the i-th column of the matrix is scaled by the i-th entry in R. This vector may be generated randomly with the "create_attack_scale_vector" function. Further, there is a "create_matrix_symmetric" function to create *symmetric* matrices (consequently with *value* zero). Finally, there is the "create_cost_vector" function for generating random row and column prices.