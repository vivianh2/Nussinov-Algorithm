import pdb
from copy import deepcopy

def nussinov_table(seq, gamma, ell=0):
    table = [[0 for _ in range(len(seq))] for _ in range(len(seq))]
    for diff in range(1+ell, len(seq)):
        for i in range(len(seq)-diff):
            j = i + diff
            opt = -1
            # case 1: if matches
            if (seq[i], seq[j]) in gamma:
                curr_opt = table[i+1][j-1] + 1
                if curr_opt > opt:
                    opt = curr_opt
            # case 2: skip head
            curr_opt = table[i+1][j]
            if curr_opt > opt:
                opt = curr_opt
            # case 3: skip tail
            curr_opt = table[i][j-1]
            if curr_opt > opt:
                opt = curr_opt
            # case 4: bifurcation
            for k in range(i+1+ell, j):
                curr_opt = table[i][k] + table[k+1][j]
                if curr_opt > opt:
                    opt = curr_opt
            table[i][j] = opt

    ans = table[0][-1]
    return ans, table

def reconstruction(seq, table, ell=0):
    n = len(seq)
    stack = [(0,n-1)]
    record = []
    while len(stack) > 0:
        i,j = stack.pop()
        print(i, j)
        if i + ell >= j:
            continue
        elif table[i+1][j] == table[i][j]:
            stack.append((i+1, j))
        elif table[i][j-1] == table[i][j]:
            stack.append((i, j-1))
        elif table[i+1][j-1] +1 == table[i][j]:
            record.append((i, j))
            stack.append((i+1, j-1))
        else:
            for k in range(i+1+ell, j):
                if table[i][k] + table[k+1][j] == table[i][j]:
                    stack.append((k+1, j))
                    stack.append((i, k))
                    break

    representation = ['.' for w in seq]
    for (i, j) in record:
        representation[i] = '('
        representation[j] = ')'
    return ''.join(representation)

def reconstruction_helper(i, j, table, ell=0):
    solutions = []
    if i + ell >= j:
        return solutions
    if table[i+1][j] == table[i][j]:
        solutions += reconstruction_helper(i+1, j, table, ell)
    if table[i][j-1] == table[i][j]:
        solutions += reconstruction_helper(i, j-1, table, ell)
    if table[i+1][j-1] + 1 == table[i][j]:
        prev_sols = reconstruction_helper(i+1, j-1, table, ell)
        if len(prev_sols) == 0:
            prev_sols = [[(i, j)]]
        else:
            prev_sols = [w+[(i,j)] for w in prev_sols]
        solutions += prev_sols
    for k in range(i+1+ell, j):
        if table[i][k] + table[k+1][j] == table[i][j]:
            left_sols = reconstruction_helper(i, k, table, ell)
            right_sols = reconstruction_helper(k+1, j, table, ell)
            combined_sols = []
            for lsol in left_sols:
                for rsol in right_sols:
                    csol = lsol + rsol
                    if len(csol) > 0:
                        combined_sols.append(csol)
            solutions += combined_sols
    return solutions

def reconstruction_all(seq, table, ell=0):
    all_sols = reconstruction_helper(0, len(seq)-1, table, ell=ell)
    reps = []
    for record in all_sols:
        representation = ['.' for w in seq]
        for (i, j) in record:
            representation[i] = '('
            representation[j] = ')'
        representation = ''.join(representation)
        reps.append(representation)
    reps = list(set(reps))
    return reps

if __name__ == "__main__":
    GAMMA = set([('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')])
    #  string = 'GCAGCAUUCG'
    #  string = 'AUAU'
    string = 'GCACUGU'
    ell = 0
    ans, table = nussinov_table(string, GAMMA, ell)
    print(ans)
    print(table)
    #  ds_format = reconstruction(string, table, ell)
    #  print(string)
    #  print(ds_format)
    reps = reconstruction_all(string, table, ell=ell)
    print('Found {} solutions'.format(len(reps)))
    for rep in reps:
        print(rep)

