# Matrix Game Designer

## Solving an instance

## Generating random instances

The "src/parameters.jl" file contains all functions for generating random matrix game instances. There are four "create_matrix" methods which all receive an entry range from which to draw entries from and a seed for reproducibility. You may pass either the number of rows and the number of columns, or create a square matrix by passing just one such quantity. Additionally, each of these has a method to pass a numerical vector R so that the i-th column of the matrix is scaled by the i-th entry in R. This vector may be generated randomly with the "create_attack_scale_vector" function. Further, there is a "create_matrix_symmetric" function to create *symmetric* matrices (consequently with *value* zero). Finally, there is the "create_cost_vector" function for generating random row and column prices.