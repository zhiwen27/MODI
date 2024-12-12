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
nbox = 40
# initial strethes in x and y direction
rho = 0.5
# final streches in x and y direction with respect to the initial
p =  0.25
# the angle the wavefunction rotates
theta = np.pi / 12
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
    cnt_supply = 0 # initialize supply counting
    cnt_demand = 0 # initialize demand counting
    for i in supply:
        if i == 0: # if the supply cell gets to zero
            cnt_supply += 1
    for i in demand:
        if i == 0: # if the demand cell gets to zero
            cnt_demand += 1
    return cnt_supply + cnt_demand # return the total number of cells

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

    return solution, cost # return solution matrix, cost

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