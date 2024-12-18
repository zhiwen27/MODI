Systems exhibiting sensitive dependence on initial conditions are defined to be chaotic. To find quantum chaos, we need to calculate the distance between quantal states and check if the difference grows exponentially. To define the distance between quantal states, we used the Wigner formulation where the states are represented in normalized real valued functions. By treating the wave functions as boxes of sands, we then need to find the minimum cost to move the sands in initial state to final state. This becomes a transportation problem that MODI method could work. <br />

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
7. If not, repeat the procedure until an optimal solution is found.

Results: <br />
![image](https://github.com/user-attachments/assets/bdaafa23-299c-425a-bd5d-312d63ae1e41) <br />
The code is able to provide answers that were close to the analytical result within 0.3% ~ 2.2% of error, and get to nbox values like 50.

Future improvements:
1. Normalizing sum_supply and sum_demand would cause a discrepency when sum_supply differs from sum_demand greatly.
2. In degeneracy the cost cells may have the same value, so I could be deleting the first value in sort list but repeatedly checking the cells.
3. Code hangs at: <br />
   nbox = 24/ 28/ 30, theta = 𝜋/3 (2 small amount of leftovers in demand[] and supply[]); <br />
   nbox = 40, theta = 𝜋/7 (print repeated values at some point); <br />
   nbox = 40, theta = 𝜋/11 (caught in an infinite loop at some point); <br />
