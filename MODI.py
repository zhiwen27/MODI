import numpy as np
import scipy as sp

# ❗️❗️ can read from a file
print("Input the rows of the cost matrix:")
row = int(input())
print("Input the columns of the cost matrix:")
col = int(input())

costmatrix = [[-1 for x in range(col)] for y in range(row)] # initialize costmatrix
costmatrix_copy = [[-1 for x in range(col)] for y in range(row)] # initialize costmatrix_copy
supply = [] # initialize supply list
demand = [] # initialize demand list
solution = [[0 for x in range(col)] for y in range(row)] # initialize solution matrix

print("Input each entry of the cost matrix:")
for i in range(0,row):
    for j in range(0,col):
        a = int(input())
        costmatrix[i][j] = a
        costmatrix_copy[i][j] = a

print("Input each entry of supply:")
for i in range(0,row):
    a = int(input())
    supply.append(a)
print("Input each entry of demand:")
for i in range(0,col):
    a = int(input())
    demand.append(a)

# using least cost method, check if all the supply and demand cells get to 0
def finished_least_cost(list):
    check = False
    for i in list:
        if i != 0: # if there's a cell in supply/ demand that is not zero, continue searching
            check = True
            break
    return check

# get u's and v's, check if all the u's and v's get their values
def finish_loop(list):
    flag = False
    for i in list:
        if i == -1e9: # having u's and v's with their initialized value
            flag = True
    return flag

# print solution in a matrix format
def print_solution(solution):
    for i in range(0,row):
        for j in range(0,col):
            print(solution[i][j],end = " ")
        print()

# print solution matrix along with u's and v's
def print_sol_u_v(solution,u_s,v_s):
    print("The solution:")
    print_solution(solution)
    print("Finding u's and v's:")
    for i in range(0,len(u_s)):
        print("u" + str(i + 1) + ": " + str(u_s[i]),end = " ")
    print()
    for i in range(0,len(v_s)):
        print("v" + str(i + 1) + ": " + str(v_s[i]),end = " ")
    print()

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
        if finished_least_cost(supply) == False and finished_least_cost(demand) == False: # if all the cells in supply and demand get to 0, end the loop
            check = False

    check_solution(costmatrix,solution,cost) # check if the solution is optimal

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
    opportunity_cost = [[1e9 for x in range(col)] for y in range(row)] # initialize opportunity cost
    solution_epsilon = [[False for x in range(col)] for y in range(row)] # marking basic cells
    u_s = [-1e9 for x in range(row)] # initialize u's
    u_s[0] = 0 # set u[0] = 0 for convenience
    v_s = [-1e9 for x in range(col)] # initialize v's
    epsilon = 1.e-9 # assign epsilon value
    continue_loop = True

    cnt_sol = check_basic_variable_sol(solution) # counting basic cells
    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0:
                solution_epsilon[i][j] = True # marking basic cells

    if cnt_sol + 1 < row + col: # if degenaracy occurs
        sort_list = []
        for i in range(0,row):
            for j in range(0,col):
                if solution[i][j] == 0: # appending non-basic cost matrix values
                    sort_list.append(costmatrix[i][j])

        sort_list.sort() # sort those non-basic cost matrix values

        while cnt_sol + 1 < row + col: # while degenaracy occurs
            for i in range(0,row):
                cnt_sol = check_basic_variable_sol(solution)
                if cnt_sol + 1 == row + col: # if the number of basic cells satisfies row + col - 1, break
                    break
                for j in range(0,col):
                    used_rows_loop = [] # initialize used_rows
                    used_cols_loop = [] # initialize used_cols
                    target_col_loop = j # mark the target column
                    # start from the smallest value and checking for non-basic cells: if no loop is formed, mark this cell to be an additional basic cell (❗️ repeated values ?)
                    if costmatrix[i][j] == sort_list[0] and solution_epsilon[i][j] == False and check_loop(solution,used_rows_loop,used_cols_loop,i,j,target_col_loop) == False:
                        solution_epsilon[i][j] = True # marking a basic cell
                        solution[i][j] += epsilon # add epsilon to it, break
                        break
            sort_list.pop(0) # delete the smallest value

    while(continue_loop): # calculate u's and v's
        find_u_s(costmatrix,u_s,v_s,0,solution_epsilon) # starting from finding v (as u[0] = 0)
        if finish_loop(u_s) == False and finish_loop(v_s) == False: # if all the u's and v's are found, break
            continue_loop = False

    print_sol_u_v(solution,u_s,v_s) # print the solution with u's and v's

    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] == 0:
                opportunity_cost[i][j] = costmatrix[i][j] - u_s[i] - v_s[j] # calculate opportunity cost = C_ij - u_i - v_j

    print("The opportunity cost:") # print opportunity cost in format
    for i in range(0,row):
        for j in range(0,col):
            if opportunity_cost[i][j] == 1e9:
                print("X", end = " ")
            else:
                print(opportunity_cost[i][j], end = " ")
        print()

    solution_state = 1 # check if the solution cost is optimal

    for i in range(0,row):
        for j in range(0,col):
            if opportunity_cost[i][j] != 1e9 and opportunity_cost[i][j] < 0: # if there's opportunity cost with negative value:
                solution_state = -1                                          # solution is not optimal
                break
    
    if solution_state == 1: # if the solution is optimal: print the cost
        print("The cost is: " + str(cost))
        print("This is the optimal solution.")

    elif solution_state == -1:
        print("The cost is: " + str(cost))
        print("This is not the optimal solution.")
        oppor_cost_min = np.min(opportunity_cost) # pivot from the lowest negative opportunity cost
        for i in range(0,row):
            for j in range(0,col):
                if opportunity_cost[i][j] == oppor_cost_min:
                    find_loop(costmatrix,i,j,solution) # find a closed loop
                    return

# calculate the number of basic variables: when solution cell is not 0
def check_basic_variable_sol(solution):
    cnt_basic_var = 0
    for i in range(0,row):
        for j in range(0,col):
            if solution[i][j] != 0: # if solution is not 0, count the cell as a basic variable
                cnt_basic_var += 1
    return cnt_basic_var # return the number of basic variables

# check if an element is already in a list
def check_list(a,list):
    check = True
    for i in list:
        if a == i:
            check = False
            break
    return check

# check if there's a loop, if not, assign the cell to be a new basic variable
# (as only checking if a loop can be formed, don't need plus_nodes and minus_nodes here)
def check_loop(solution,used_rows,used_cols,start_row,start_col,target_col):
    loop_found = False # set loop_found to false
    check_loop_row(solution,used_rows,used_cols,start_row,start_col,target_col,loop_found) # check if there's a loop starting from column
    return loop_found # return if a loop is found

# search along a row and check if there's basic cell that could form a loop
def check_loop_row(solution,used_rows,used_cols,start_row,start_col,target_col,loop_found):
    end_loop = False # check if reaching the end of the loop
    len_used_rows = len(used_rows)
    for i in range(0,col): # searching along all the columns on the row
        if loop_found == True: # if a loop is found, return
            return
        # if there's another basic cell (excluding start_col itself) on the row, and this row hasn't been visited before
        if i != start_col and solution[start_row][i] != 0 and check_list(start_row,used_rows):
            used_rows.append(start_row) # add this row to used_rows
            if target_col == i: # if the column now equals to the target_col, a loop forms, set loop_found to true and return
                loop_found = True
                return
            else: # otherwise, starting from the new cell and search along the column
                check_loop_col(solution,used_rows,used_cols,start_row,i,target_col,loop_found)
        if i == (col - 1) and len(used_rows) == len_used_rows: # if reaching the end of the loop and no new elements have been added to used_rows, then check_loop_row fails
            end_loop = True # having reached the end of the loop
    if end_loop == True: # if check_loop_row fails, return
        if loop_found == True: # if a loop is formed by the last element in the loop (there's no chance looking back in the loop), return
            return
        if len(used_cols) != 0: # delete the last element in used_cols
            used_cols.pop()
        return

# search along a column and check if there's basic cell that could form a loop
def check_loop_col(solution,used_rows,used_cols,start_row,start_col,target_col,loop_found):
    end_loop = False # check if reaching the end of the loop
    len_used_cols = len(used_cols)
    for i in range(0,row): # searching all the rows on the column
        if loop_found == True: # if a loop is found, return
            return
        # if there's another basic cell (excluding start_row itself) on the column, and this column hasn't been visited before
        if i != start_row and solution[i][start_col] != 0 and check_list(start_col,used_cols):
            used_cols.append(start_col) # add this column to used_cols
            check_loop_row(solution,used_rows,used_cols,i,start_col,target_col,loop_found) # starting from the new cell and search along the row
        if i == (row - 1) and len(used_cols) == len_used_cols: # if reaching the end of the loop and no new elements have been added to used_cols, then check_loop_col fails
            end_loop = True
    if end_loop == True: # if check_loop_col fails, return
        if loop_found == True:
            return
        if len(used_rows) != 0: # delete the last element in used_rows
            used_rows.pop()
        return

# search along a row and find if there's basic cell to be the next node on the loop
def search_row(solution,plus_nodes,minus_nodes,used_rows,used_cols,start_row,start_col,target_col,loop_found):
    end_loop = False # check if reaching the end of the loop
    len_minus = len(minus_nodes)
    for i in range(0,col): # searching along all the columns on the row
        if loop_found == True: # if a loop is found, return
            return
        # if there's another basic cell (excluding start_col itself) on the row, and this row hasn't been visited before
        if i != start_col and solution[start_row][i] != 0 and check_list(start_row,used_rows):
            minus_nodes.append([start_row,i]) # add this cell to minus_nodes
            used_rows.append(start_row) # add this row to used_rows
            if target_col == i: # if the column now equals to the target_col, a loop forms, set loop_found to true and return
                loop_found = True
                return
            else: # otherwise, starting from the new cell and search along the column
                search_col(solution,plus_nodes,minus_nodes,used_rows,used_cols,start_row,i,target_col,loop_found)
        if i == (col - 1) and len(minus_nodes) == len_minus: # if reaching the end of the loop and no new elements have been added to minus_nodes, then search_col fails
            end_loop = True
    if end_loop == True: # if search_row fails, return
        if loop_found == True: # if a loop is formed by the last element in the loop, return
            return
        plus_nodes.pop() # delete the last element in plus_nodes
        used_cols.pop() # delete the last element in used_cols
        return

# search along a column and find if there's basic cell to be the next node on the loop
def search_col(solution,plus_nodes,minus_nodes,used_rows,used_cols,start_row,start_col,target_col,loop_found):
    end_loop = False # check if reaching the end of the loop
    len_plus = len(plus_nodes)
    for i in range(0,row): # searching all the rows on the column
        if loop_found == True: # if a loop is found, return
            return
        # if there's another basic cell (excluding start_row itself) on the column, and this column hasn't been visited before
        if i != start_row and solution[i][start_col] != 0 and check_list(start_col,used_cols):
            plus_nodes.append([i,start_col]) # add this cell to plus_nodes
            used_cols.append(start_col) # add this col to used_cols
            search_row(solution,plus_nodes,minus_nodes,used_rows,used_cols,i,start_col,target_col,loop_found) # starting from the new cell and search along the row
        if i == (row - 1) and len(plus_nodes) == len_plus: # if reaching the end of the loop and no new elements have been added to plus_nodes, then search_row fails
            end_loop = True
    if end_loop == True: # if search_col fails, return
        if loop_found == True:
            return
        minus_nodes.pop() # delete the last element in minus_nodes
        used_rows.pop() # delete the last element in used_rows
        return

# find a closed loop: start with the unoccupied cell with largest negative opportunity cost and mark this cell with a plus sign (1),
# trace a path along the row to an occupied cell (solution != 0), and mark with a minus sign (-1), then continue down the column and
# repeat the process till it goes back to the starting cell
def find_loop(costmatrix,start_row,start_col,solution):
    closed_loop = [[0 for x in range(col)] for y in range(row)] # initialize closed_loop
    closed_loop[start_row][start_col] = 1 # mark the cell with the lowest opportunity cost to be the starting point

    plus_nodes = [] # initialize plus_nodes: store the position of each cell on the loop
    plus_nodes.append([start_row,start_col]) # add the pivot cell to plus_nodes list
    minus_nodes = [] # initialize minus_nodes: store the position of each cell on the loop
    used_rows = [] # initialize used_rows
    used_cols = [] # initialize used_cols
    target_col = start_col # mark the target column

    loop_found = False # set loop_found to false first

    # starting from the pivot cell and search along the row
    search_row(solution,plus_nodes,minus_nodes,used_rows,used_cols,start_row,start_col,target_col,loop_found)
    
    for a in plus_nodes:
        closed_loop[a[0]][a[1]] = 1 # mark all the cells in plus_nodes to be 1

    for a in minus_nodes:
        closed_loop[a[0]][a[1]] = -1 # mark all the cells in minus_nodes to be -1

    print("Marking a closed loop:") # print the closed loop in a matrix format
    for i in range(0,row):
        for j in range(0,col):
            print(closed_loop[i][j],end = " ")
        print()

    cost = 0 # set cost back to 0
    min_sandmove = 1e9 # initialize min_sandmove with a max value
    epsilon = 1.e-9 # set epsilon to 1.e-9

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
            if solution[i][j] != epsilon: # calculate total cost: excluding the cells equal to epsilon (❗️ not accurate: the cells plus/ minus epsilon are included)
                cost += costmatrix[i][j] * solution[i][j]
    
    check_solution(costmatrix,solution,cost) # check if solution is optimal

find_initial_sol(costmatrix,costmatrix_copy,supply,demand) # start from finding initial solution
