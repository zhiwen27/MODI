This Python code uses MODI Method to solve Transportation Problem occurred in Quantum Chaos:
MODI Method:
1. Find an initial solution with m + n -1 occupied cells and calculate u[i] and v[j] for rows and columns using Least Cost Method.
To start with, assign u[i] to be zero. Then calculate of u_s and v_s using C_ij = u_i + v_j. 
2. For unoccupied cells, calculate opportunity cost = C_ij - u_i - v_j. 
3. Test if basic feasible solution is optimal:
   (1). Check whether row(m) + col(n) - 1 = total number of allocated cells. If the total number of allocated cells is fewer than m + n - 1, the case of degeneracy occurs. To convert unoccupied cells into occupied cells, start from the least value among all the unoccupied cells and check if that cell has no closed-loop formation. If so, select it as a new occupied cell and assign with value epsilon.
   (2). Examine the sign of each opportunity cost cell:
        If opportunity_cost > 0, then current solution is optimal.
        If opportunity_cost = 0, then an alternative solution exists; If there exists opportunity_cost < 0, then an improved solution can be obtained by pivoting from the cell with the largest negative value of opportunity cost.
4. Construct a closed loop from the cell with largest negative opportunity cost. Start the loop with the selected cell and mark it
with a plus sign. Then, search along the row to an occupied cell, mark the cell with a minus sign and continue down the column to an
occupied cell and mark the corner with a plus sign. Continue the search until the loop gets back to the selected cell.
5. Find the smallest value among the cells marked with minus sign. Add it to cells marked with plus sign and subtract it from cells
marked with minus sign.
6. Check again for the revised solution to see if it's optimal.
Quantum Chaos:


Errors:
1. Normalizing sum_supply and sum_demand would cause a discrepency when sum_supply differs from sum_demand greatly.
2. Codes hang at nbox = 24/ 28/ 30, theta = 𝜋/3 (2 small amount of leftouts in demand[] and supply[]);
                 nbox = 40, theta = 𝜋/7 (print repeated values at some point);
                 nbox = 40, theta = 𝜋/11 (hang at some point);

References:
