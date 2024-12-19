import sys
sys.setrecursionlimit(5000)  # set maximum recursion depth manually

import numpy as np
import scipy as sp
import time as tm
from scipy.spatial import distance

# set parameters of wavefunctions to be defined
# xmin, xmax, ymin, ymax set the border of the sandbox
xmin = -3.747
xmax = 3.747
ymin = -3.747
ymax = 3.747
# nbox sets the number of grids on one side. The sandbox will be divided into nbox * nbox grids.
nbox = 25
# initial strethes in x and y direction
rho = 0.5
# final streches in x and y direction with respect to the initial
p =  0.25
# the angle the wavefunction rotates
theta = np.pi / 11
# set cutoff value to be 0.001 / (nbox^2) (arbitrarily)
gridboxcutoff = 0.001 / (nbox ** 2)

# define initial and final wavefunctions
def f_init(y,x):
    return (1/np.pi)*np.e**((-x**2)-((y**2)/(rho**2)))

aa = ((np.cos(theta)**2)/p**2)+(((p**2)*(np.sin(theta)**2))/(rho**2))
bb = 2*np.sin(theta)*np.cos(theta)*((1/p**2)-(p**2/rho**2))
cc = (np.sin(theta)**2/p**2)+(p**2*np.cos(theta)**2/rho**2)

def f_final(y,x):
    return (1/np.pi)*np.e**(-((aa*(x**2))+(bb*x*y)+(cc*(y**2))))

init_integral = sp.integrate.dblquad(f_init, xmin, xmax, lambda x: ymin, lambda x: ymax)
fin_integral = sp.integrate.dblquad(f_final, xmin, xmax, lambda x: ymin, lambda x: ymax)

# normalize initial and final wavefunctions
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

# normalize sum_supply and sum_demand
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
        if i == -1e9: # having u's and v's with their initialized value
            flag = True
    return flag

# find the initial solution using least cost method
def find_initial_sol(costmatrix,costmatrix_copy,supply,demand):
    cost = 0 # initialize cost
    matrix_min = np.min(costmatrix_copy) # starting from the minimum
    max_val = 1e9 # define a max value

    check = True

    while(check):
        for i in range(0,row):
            for j in range(0,col):
                if costmatrix_copy[i][j] == matrix_min: # if costmatrix_copy equals the minimum
                    if supply[i] >= demand[j]: # if supply is bigger than demand
                        supply[i] -= demand[j] # subtract the same demand value from corresponding supply cell
                        solution[i][j] += demand[j] # add the demand value to the solution cell
                        cost += costmatrix[i][j] * demand[j] # add this step to the total cost
                        demand[j] = 0 # remove all the demand value
                    else: # if demand is bigger than supply
                        demand[j] -= supply[i] # subtract the same supply value from corresponding demand cell
                        solution[i][j] += supply[i] # add the supply value to the solution cell
                        cost += costmatrix[i][j] * supply[i] # add this step to the total cost
                        supply[i] = 0 # remove all the supply value
                    costmatrix_copy[i][j] = max_val # label this costmatrix_copy cell with the defined max value
                    matrix_min = np.min(costmatrix_copy) # check the minimum in the rest cells
        # if all the cells in supply and demand get to 0 or there's one cell remaining some very small amount of sand, end the loop
        if finished_least_cost(supply,demand) == len(supply) + len(demand) - 1 or finished_least_cost(supply,demand) == len(supply) + len(demand):
            check = False

    return solution, cost

# find v's, start from row 0
def find_u_s(costmatrix,u_s,v_s,sol_row,solution_epsilon):
    for j in range(0,col):
        if solution_epsilon[sol_row][j] == True: # for a corresponding basic cell, calculate v value (start from u[0] = 0)
            if v_s[j] == -1e9:
                v_s[j] = costmatrix[sol_row][j] - u_s[sol_row]
                find_v_s(costmatrix,u_s,v_s,j,solution_epsilon) # continue searching for u value

# find u's
def find_v_s(costmatrix,u_s,v_s,sol_col,solution_epsilon):
    for i in range(0,row):
        if solution_epsilon[i][sol_col] == True: # for a corresponding basic cell, calculate u value
            if u_s[i] == -1e9:
                u_s[i] = costmatrix[i][sol_col] - v_s[sol_col]
                find_u_s(costmatrix,u_s,v_s,i,solution_epsilon) # continue searching for v value

# assign basic variables, find u's and v's, calculate opportunity cost and decide if the cost is the optimal solution
def check_solution(costmatrix,solution,cost):
    opportunity_cost = [[1e9 for x in range(col)] for y in range(row)]
    solution_epsilon = [[False for x in range(col)] for y in range(row)] # marking basic cells
    u_s = [-1e9 for x in range(row)]
    u_s[0] = 0 # set u[0] = 0 for convenience
    v_s = [-1e9 for x in range(col)]
    epsilon = 1.e-9 # assign epsilon value
    continue_loop = True

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0:
                solution_epsilon[i][j] = True # marking basic cells

    cnt_sol = check_basic_variable_sol(solution) # counting basic cells

    if cnt_sol + 1 < row + col: # if degenaracy occurs
        sort_list = []
        for i in range(0,row):
            for j in range(0,col):
                if solution[i][j] == 0: # appending non-basic cost matrix values
                    sort_list.append(costmatrix[i][j])

        sort_list.sort() # sort the non-basic cost matrix values
        while cnt_sol + 1 < row + col: # while degenaracy occurs
            for i in range(0,row):
                cnt_sol = check_basic_variable_sol(solution)
                if cnt_sol + 1 == row + col: # if the number of basic cells satisfies row + col - 1, break
                    break
                for j in range(0,col):
                    used_rows_loop = []
                    used_cols_loop = []
                    target_col_loop = j # mark the target column
                    # start from the smallest value and checking for non-basic cells: if no loop is formed, mark this cell to be an additional basic cell
                    if costmatrix[i][j] == sort_list[0] and solution_epsilon[i][j] == False and check_loop(solution,used_rows_loop,used_cols_loop,i,j,target_col_loop) == False:
                        solution_epsilon[i][j] = True # marking a basic cell
                        solution[i][j] += epsilon # add epsilon to it, break
                        break
            sort_list.pop(0) # delete the smallest value
    
    while(continue_loop): # calculate u's and v's
        find_u_s(costmatrix,u_s,v_s,0,solution_epsilon) # starting from finding v (as u[0] = 0)
        if finish_loop(u_s) == False and finish_loop(v_s) == False: # if all the u's and v's are found, break
            continue_loop = False

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] == 0:
                opportunity_cost[i][j] = costmatrix[i][j] - u_s[i] - v_s[j] # calculate opportunity cost = C_ij - u_i - v_j

    solution_state = 1 # check if the solution cost is optimal

    for i in range(0,row):
        for j in range(0,col):
            if opportunity_cost[i][j] != 1e9 and opportunity_cost[i][j] < 0: # if there's opportunity cost with negative value:
                solution_state = -1                                          # solution is not optimal
                break
    
    if solution_state == 1: # if the solution is optimal: print the cost
        print("The cost is: " + str(cost))
        print("This is the optimal solution.")
        return False, solution, cost

    elif solution_state == -1:
        print("The cost is: " + str(cost))
        print("This is not the optimal solution.")
        oppor_cost_min = np.min(opportunity_cost) # pivot from the lowest negative opportunity cost
        for i in range(0,row):
            for j in range(0,col):
                if opportunity_cost[i][j] == oppor_cost_min:
                    solution, cost = find_loop(costmatrix,i,j,solution) # find a closed loop
                    return True, solution, cost

# calculate the number of basic variables: when solution cell is not 0
def check_basic_variable_sol(solution):
    cnt_basic_var = 0
    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0: # if solution is not 0, count the cell as a basic variable
                cnt_basic_var += 1
    return cnt_basic_var

# check if an element is already in a list
def check_list(a,list):
    check = True
    for i in list:
        if a == i:
            check = False
            break
    return check

# check if there's a loop, if not, assign the cell to be a new basic variable (together with check_loop_run)
def check_loop(solution,used_rows,used_cols,start_row,start_col,pivotcol):
    loop_found = False
    check_loop_run(solution,used_rows,used_cols,start_row,start_col,loop_found,pivotcol)
    return loop_found

# generally follows the logic of find_loop
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

# find a closed loop: start with the unoccupied cell with largest negative opportunity cost and mark this cell with a plus sign (1),
# trace a path along the row to an occupied cell (solution != 0), and mark with a minus sign (-1), then continue down the column and
# repeat the process till it goes back to the starting cell
def find_loop(costmatrix,start_row,start_col,solution):
    closed_loop = [[0 for x in range(col)] for y in range(row)]
    closed_loop[start_row][start_col] = 1 # mark the cell with the lowest opportunity cost to be the starting point

    plus_nodes = []
    minus_nodes = []
    used_rows = []
    used_cols = []
    pivotcol = start_col # mark the target column
    loop = [[start_row,start_col]] # add the pivot cell to plus_nodes list
    loop_found = False
    hor = 0 # set horizontal search flag
    vert = 1 # set vertical search flag
    search = hor # start by searching horizontally
    currentrow = start_row # row of current node (initialize with start_row)
    currentcol = start_col # column of current node (initialize with start_column)
    startd = 0 # index of where to start horizontal search (usually starts at 0)
    starts = 0 # index of where to start vertical search (usually starts at 0)

    while(loop_found == False):
        if search == hor: # horizontal search
            for d in range(startd,col):
                # if there's another basic cell (excluding start_col itself) on the row, and this row hasn't been visited before
                if d != currentcol and solution[currentrow][d] != 0 and check_list(d,used_cols):
                    loop.append([currentrow,d]) # add this cell to loop_nodes
                    used_rows.append(currentrow) # add this row to used_rows
                    if d == pivotcol: # if the column now equals to the target_col, a loop forms, set loop_found to true
                        loop_found = True
                    search = vert # switch to vertical search
                    currentcol = d # mark current column
                    break
            startd = 0 # set the index starting horizontal search back to 0
            if search == hor: # if finishes searching the row without finding a usable node, resume the last vertical search
                if len(loop) == 1: # no usable node on the initial row, the search fails
                    break
                starts = loop[len(loop) - 1][0] + 1 # start from the next row index
                loop.pop() # delete the last element in loop
                used_cols.pop() # delete the last element in used_cols
                currentrow = loop[len(loop) - 1][0] # mark current row
                currentcol = loop[len(loop) - 1][1] # mark current column
                search = vert # switch to vertical search
        if loop_found == False and search == vert:
            for s in range(starts,row):
                # if there's another basic cell (excluding start_row itself) on the column, and this column hasn't been visited before
                if s != currentrow and solution[s][currentcol] != 0 and check_list(s,used_rows):
                    loop.append([s,currentcol]) # add this cell to loop_nodes
                    used_cols.append(currentcol) # add this column to used_cols
                    search = hor # switch to horizontal search
                    currentrow = s # mark current row
                    break
            starts = 0 # set the index starting vertical search back to 0
            if search == vert: # if finishes searching the column without finding a usable node, resume the last horizontal search
                startd = loop[len(loop) - 1][1] + 1 # start from the next column index
                loop.pop() # delete the last element in loop
                used_rows.pop() # delete the last element in used_rows
                currentrow = loop[len(loop) - 1][0] # mark current row
                currentcol = loop[len(loop) - 1][1] # mark current column
                search = hor # switch to horizontal search
    
    # separate loop nodes to plus nodes and minus nodes
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
    min_sandmove = 1e9 # initialize min_sandmove with a MAX value
    epsilon = 1.e-9 # assign epsilon value

    for i in range(0,row):
        for j in range(0,col):
            if closed_loop[i][j] == -1: # find the smallest value among minus_nodes and set min_sandmove to corresponding solution value
                if solution[i][j] < min_sandmove:
                    min_sandmove = solution[i][j]
    
    for i in range(0,row):
        for j in range(0,col):
            if closed_loop[i][j] == -1: # subtract min_sandmove from the occupied cells marked with minus sign
                solution[i][j] -= min_sandmove
            elif closed_loop[i][j] == 1: # add min_sandmove to the occupied cells marked with plus sign
                solution[i][j] += min_sandmove

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != epsilon:  # calculate total cost (excluding the cells equal to epsilon for better precision)
                cost += costmatrix[i][j] * solution[i][j]
    
    return solution, cost

start_time = tm.time() # mark start time
solution, cost = find_initial_sol(costmatrix,costmatrix_copy,supply,demand) # find an initial solution
continue_check_sol, solution, cost = check_solution(costmatrix,solution,cost)
while(continue_check_sol): # keep on finding a better solution until the optimal is found
    continue_check_sol, solution, cost = check_solution(costmatrix,solution,cost)

end_time = tm.time() # mark end time
print(end_time - start_time) # print out needed time
