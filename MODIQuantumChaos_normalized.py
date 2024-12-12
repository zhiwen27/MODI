'''
Using MODI method to solve transportation problem:
1. Find an initial solution with m + n -1 occupied cells and calculate u[i] and v[j] for rows and columns using Least Cost Method.
To start with, assign u[i] to be zero. Then calculate of u_s and v_s using C_ij = u_i + v_j. 
2. For unoccupied cells, calculate opportunity cost = C_ij - u_i - v_j. 
3. Test if basic feasible solution is optimal: 1). Check whether row(m) + col(n) - 1 = total number of allocated cells. If the total
number of allocated cells is fewer than m + n - 1, the case of degeneracy occurs. To convert unoccupied cells into occupied cells,
start from the least value among all the unoccupied cells and check if that cell has no closed-loop formation. If so, select it as a
new occupied cell and assign with value epsilon. 2). Examine the sign of each opportunity cost cell: If opportunity_cost > 0, then
current solution is optimal; If opportunity_cost = 0, then an alternative solution exists; If there exists opportunity_cost < 0, then
an improved solution can be obtained by pivoting from the cell with the largest negative value of opportunity cost.
4. Construct a closed loop from the cell with largest negative opportunity cost. Start the loop with the selected cell and mark it
with a plus sign. Then, search along the row to an occupied cell, mark the cell with a minus sign and continue down the column to an
occupied cell and mark the corner with a plus sign. Continue the search until the loop gets back to the selected cell.
5. Find the smallest value among the cells marked with minus sign. Add it to cells marked with plus sign and subtract it from cells
marked with minus sign.
6.Check again for the revised solution to see if it's optimal.

Errors:
1. The normalizing method would cause a discrepency when sum_supply differs from sum_demand greatly.
2. Hang at nbox = 24/ 28/ 30, theta = np.pi / 3 (2 small amount of leftouts in demand[] and supply[]).
'''
import sys
sys.setrecursionlimit(5000)  # set maximum recursion depth manually

import numpy as np
import scipy as sp
import time as tm
from scipy.spatial import distance

# input parameters
xmin = -3.747
xmax = 3.747
ymin = -3.747
ymax = 3.747
nbox = 50
rho = 0.5
p =  0.25
theta = np.pi / 3
gridboxcutoff = 0.001 / (nbox ** 2)

# define initial and final functions
def f_init(y,x):
    return (1/np.pi)*np.e**((-x**2)-((y**2)/(rho**2)))

aa = ((np.cos(theta)**2)/p**2)+(((p**2)*(np.sin(theta)**2))/(rho**2))
bb = 2*np.sin(theta)*np.cos(theta)*((1/p**2)-(p**2/rho**2))
cc = (np.sin(theta)**2/p**2)+(p**2*np.cos(theta)**2/rho**2)

def f_final(y,x):
    return (1/np.pi)*np.e**(-((aa*(x**2))+(bb*x*y)+(cc*(y**2))))

init_integral = sp.integrate.dblquad(f_init, xmin, xmax, lambda x: ymin, lambda x: ymax)
fin_integral = sp.integrate.dblquad(f_final, xmin, xmax, lambda x: ymin, lambda x: ymax)

def f_init_norm(y,x):
    return f_init(y,x)/init_integral[0]

def f_final_norm(y,x):
    return f_final(y,x)/fin_integral[0]

def integrand(y,x):
    return f_init_norm(y,x)-f_final_norm(y,x)

# calculate sand boxes, supply and demand cells
sand = [[0 for x in range(nbox)] for y in range(nbox)]
dx = (xmax-xmin)/nbox
dy = (ymax-ymin)/nbox
for i in range(nbox):
    for j in range(nbox):
        diff = sp.integrate.dblquad(integrand,xmin+i*dx,xmin+(i+1)*dx,lambda x: ymin+j*dy,lambda x: ymin+(j+1)*dy)
        sand[i][j] = diff[0]

outboxes = []
supply = []
inboxes = []
demand_store = []
demand = []
for i in range(nbox):
    for j in range(nbox):
        if sand[i][j] > gridboxcutoff:
            outboxes.append([xmin+(i-.5)*dx,ymin + (j - .5)*dy])
            supply.append(sand[i][j])
        elif sand[i][j] < -gridboxcutoff:
            inboxes.append([xmin+(i-.5)*dx,ymin + (j - .5)*dy])
            demand_store.append(sand[i][j])
demand = [abs(i) for i in demand_store]

sum_demand = 0
for i in demand:
    sum_demand += i
sum_supply = 0
for i in supply:
    sum_supply += i

for i in range(0,len(supply)):
    supply[i] = supply[i] * (sum_supply + sum_demand) / (2 * sum_supply)
for i in range(0,len(demand)):
    demand[i] = demand[i] * (sum_supply + sum_demand) / (2 * sum_demand)

# find cost matrices
row = len(supply)
col = len(demand)
costmatrix = [[-1 for x in range(col)] for y in range(row)]
costmatrix_copy = [[-1 for x in range(col)] for y in range(row)]
solution = [[0 for x in range(col)] for y in range(row)]
for i in range(len(supply)):
    for j in range(len(demand)):
        euc_dist = distance.euclidean(inboxes[j],outboxes[i])
        costmatrix[i][j] = euc_dist
        costmatrix_copy[i][j] = euc_dist

# using least cost method, count the total number of supply and demand cells that get to 0
def finished_least_cost(supply,demand):
    cnt_supply = 0
    cnt_demand = 0
    for i in supply:
        if i == 0:
            cnt_supply += 1
    for i in demand:
        if i == 0:
            cnt_demand += 1
    return cnt_supply + cnt_demand

# get u's and v's, check if all the u's and v's get their values
# input a list, return true if there's still u/ v with its initialized value
def finish_loop(list):
    flag = False
    for i in list:
        if i == -1e9:
            flag = True
    return flag

# find the initial solution using least cost method
def find_initial_sol(costmatrix,costmatrix_copy,supply,demand):
    cost = 0
    matrix_min = np.min(costmatrix_copy)
    max_val = 1e9

    check = True

    while(check):
        for i in range(0,row):
            for j in range(0,col):
                if costmatrix_copy[i][j] == matrix_min:
                    if supply[i] >= demand[j]:
                        supply[i] -= demand[j]
                        solution[i][j] += demand[j]
                        cost += costmatrix[i][j] * demand[j]
                        demand[j] = 0
                    else:
                        demand[j] -= supply[i]
                        solution[i][j] += supply[i]
                        cost += costmatrix[i][j] * supply[i]
                        supply[i] = 0
                    costmatrix_copy[i][j] = max_val
                    matrix_min = np.min(costmatrix_copy)
        if finished_least_cost(supply,demand) == len(supply) + len(demand) - 1 or finished_least_cost(supply,demand) == len(supply) + len(demand):
            check = False

    return solution, cost

def find_u_s(costmatrix,u_s,v_s,sol_row,solution_epsilon):
    for j in range(0,col):
        if solution_epsilon[sol_row][j] == True:
            if v_s[j] == -1e9:
                v_s[j] = costmatrix[sol_row][j] - u_s[sol_row]
                find_v_s(costmatrix,u_s,v_s,j,solution_epsilon)

def find_v_s(costmatrix,u_s,v_s,sol_col,solution_epsilon):
    for i in range(0,row):
        if solution_epsilon[i][sol_col] == True:
            if u_s[i] == -1e9:
                u_s[i] = costmatrix[i][sol_col] - v_s[sol_col]
                find_u_s(costmatrix,u_s,v_s,i,solution_epsilon)

def check_solution(costmatrix,solution,cost):
    opportunity_cost = [[1e9 for x in range(col)] for y in range(row)]
    solution_epsilon = [[False for x in range(col)] for y in range(row)]
    u_s = [-1e9 for x in range(row)]
    u_s[0] = 0
    v_s = [-1e9 for x in range(col)]
    epsilon = 1.e-9
    continue_loop = True

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0:
                solution_epsilon[i][j] = True

    cnt_sol = check_basic_variable_sol(solution)

    if cnt_sol + 1 < row + col:
        sort_list = []
        for i in range(0,row):
            for j in range(0,col):
                if solution[i][j] == 0:
                    sort_list.append(costmatrix[i][j])

        sort_list.sort()
        while cnt_sol + 1 < row + col:
            for i in range(0,row):
                cnt_sol = check_basic_variable_sol(solution)
                if cnt_sol + 1 == row + col:
                    break
                for j in range(0,col):
                    used_rows_loop = []
                    used_cols_loop = []
                    target_col_loop = j
                    if costmatrix[i][j] == sort_list[0] and solution_epsilon[i][j] == False and check_loop(solution,used_rows_loop,used_cols_loop,i,j,target_col_loop) == False:
                        solution_epsilon[i][j] = True
                        solution[i][j] += epsilon
                        break
            sort_list.pop(0)
    
    while(continue_loop):
        find_u_s(costmatrix,u_s,v_s,0,solution_epsilon)
        if finish_loop(u_s) == False and finish_loop(v_s) == False:
            continue_loop = False

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] == 0:
                opportunity_cost[i][j] = costmatrix[i][j] - u_s[i] - v_s[j]

    solution_state = 1

    for i in range(0,row):
        for j in range(0,col):
            if opportunity_cost[i][j] != 1e9 and opportunity_cost[i][j] < 0:
                solution_state = -1
                break
    
    if solution_state == 1:
        print("The cost is: " + str(cost))
        print("This is the optimal solution.")
        return False, solution, cost

    elif solution_state == -1:
        print("The cost is: " + str(cost))
        print("This is not the optimal solution.")
        oppor_cost_min = np.min(opportunity_cost)
        for i in range(0,row):
            for j in range(0,col):
                if opportunity_cost[i][j] == oppor_cost_min:
                    solution, cost = find_loop(costmatrix,i,j,solution)
                    return True, solution, cost

def check_basic_variable_sol(solution):
    cnt_basic_var = 0
    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0:
                cnt_basic_var += 1
    return cnt_basic_var

def check_list(a,list):
    check = True
    for i in list:
        if a == i:
            check = False
            break
    return check

def check_loop(solution,used_rows,used_cols,start_row,start_col,pivotcol):
    loop_found = False
    check_loop_run(solution,used_rows,used_cols,start_row,start_col,loop_found,pivotcol)
    return loop_found

def check_loop_run(solution,used_rows,used_cols,start_row,start_col,loop_found,pivotcol):
    loop = [[start_row,start_col]]
    hor = 0
    vert = 1
    search = hor
    currentrow = start_row
    currentcol = start_col
    startd = 0
    starts = 0

    while(loop_found == False):
        if search == hor:
            for d in range(startd,col):
                startd = 0
                if d != currentcol and solution[currentrow][d] != 0 and check_list(d,used_cols):
                    loop.append([currentrow,d])
                    used_rows.append(currentrow)
                    if d == pivotcol:
                        loop_found = True
                    search = vert
                    currentcol = d
                    break
            if search == hor:
                if len(loop) == 1:
                    break
                starts = loop[len(loop) - 1][0] + 1
                loop.pop()
                used_cols.pop()
                currentrow = loop[len(loop) - 1][0]
                currentcol = loop[len(loop) - 1][1]
                search = vert
        if loop_found == False and search == vert:
            for s in range(starts,row):
                starts = 0
                if s != currentrow and solution[s][currentcol] != 0 and check_list(s,used_rows):
                    loop.append([s,currentcol])
                    used_cols.append(currentcol)
                    search = hor
                    currentrow = s
                    break
            if search == vert:
                startd = loop[len(loop) - 1][1] + 1
                loop.pop()
                used_rows.pop()
                currentrow = loop[len(loop) - 1][0]
                currentcol = loop[len(loop) - 1][1]
                search = hor
    return loop_found

def find_loop(costmatrix,start_row,start_col,solution):
    closed_loop = [[0 for x in range(col)] for y in range(row)]
    closed_loop[start_row][start_col] = 1

    plus_nodes = []
    minus_nodes = []
    used_rows = []
    used_cols = []
    pivotcol = start_col
    loop = [[start_row,start_col]]
    loop_found = False
    hor = 0
    vert = 1
    search = hor
    currentrow = start_row
    currentcol = start_col
    startd = 0
    starts = 0

    while(loop_found == False):
        if search == hor:
            for d in range(startd,col):
                if d != currentcol and solution[currentrow][d] != 0 and check_list(d,used_cols):
                    loop.append([currentrow,d])
                    used_rows.append(currentrow)
                    if d == pivotcol:
                        loop_found = True
                    search = vert
                    currentcol = d
                    break
            startd = 0
            if search == hor:
                if len(loop) == 1:
                    break
                starts = loop[len(loop) - 1][0] + 1
                loop.pop()
                used_cols.pop()
                currentrow = loop[len(loop) - 1][0]
                currentcol = loop[len(loop) - 1][1]
                search = vert
        if loop_found == False and search == vert:
            for s in range(starts,row):
                if s != currentrow and solution[s][currentcol] != 0 and check_list(s,used_rows):
                    loop.append([s,currentcol])
                    used_cols.append(currentcol)
                    search = hor
                    currentrow = s
                    break
            starts = 0
            if search == vert:
                startd = loop[len(loop) - 1][1] + 1
                loop.pop()
                used_rows.pop()
                currentrow = loop[len(loop) - 1][0]
                currentcol = loop[len(loop) - 1][1]
                search = hor

    for i in range(0,len(loop)):
        if (i%2 == 0):
            plus_nodes.append(loop[i])
        else:
            minus_nodes.append(loop[i])

    for a in plus_nodes:
        closed_loop[a[0]][a[1]] = 1

    for a in minus_nodes:
        closed_loop[a[0]][a[1]] = -1

    cost = 0
    min_sandmove = 1e9
    epsilon = 1.e-9

    for i in range(0,row):
        for j in range(0,col):
            if closed_loop[i][j] == -1:
                if solution[i][j] < min_sandmove:
                    min_sandmove = solution[i][j]
    
    for i in range(0,row):
        for j in range(0,col):
            if closed_loop[i][j] == -1:
                solution[i][j] -= min_sandmove
            elif closed_loop[i][j] == 1:
                solution[i][j] += min_sandmove

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != epsilon:
                cost += costmatrix[i][j] * solution[i][j]
    
    return solution, cost

 # mark start time
start_time = tm.time()
solution, cost = find_initial_sol(costmatrix,costmatrix_copy,supply,demand)
continue_check_sol, solution, cost = check_solution(costmatrix,solution,cost)
while(continue_check_sol):
    continue_check_sol, solution, cost = check_solution(costmatrix,solution,cost)

end_time = tm.time()
print(end_time - start_time)