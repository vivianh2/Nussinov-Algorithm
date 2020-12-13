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

def reconstruction(seq, table, gamma, ell=0):
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

    return record

    #  representation = ['.' for w in seq]
    #  for (i, j) in record:
    #      representation[i] = '('
    #      representation[j] = ')'
    #  return ''.join(representation)

def reconstruction_helper(seq, i, j, table, gamma, ell=0):
    if i == 0 and j == 2:
        pdb.set_trace()
    solutions = []
    if i + ell >= j:
        return solutions
    if table[i+1][j] == table[i][j]:
        solutions += reconstruction_helper(seq, i+1, j, table, gamma, ell)
    if table[i][j-1] == table[i][j]:
        solutions += reconstruction_helper(seq, i, j-1, table, gamma, ell)
    if table[i+1][j-1] + 1 == table[i][j] and (seq[i], seq[j]) in gamma:
        prev_sols = reconstruction_helper(seq, i+1, j-1, table, gamma, ell)
        if len(prev_sols) == 0:
            prev_sols = [[(i, j)]]
        else:
            prev_sols = [w+[(i,j)] for w in prev_sols]
        solutions += prev_sols
    for k in range(i+1+ell, j):
        if table[i][k] + table[k+1][j] == table[i][j]:
            left_sols = reconstruction_helper(seq, i, k, table, gamma, ell)
            right_sols = reconstruction_helper(seq, k+1, j, table, gamma, ell)
            combined_sols = []
            for lsol in left_sols:
                for rsol in right_sols:
                    csol = lsol + rsol
                    if len(csol) > 0:
                        combined_sols.append(csol)
            solutions += combined_sols
    return solutions

def reconstruction_all(seq, table, gamma, ell=0):
    all_sols = reconstruction_helper(seq, 0, len(seq)-1, table, gamma, ell=ell)
    for i in range(len(all_sols)):
        record = all_sols[i]
        all_sols[i] = sorted(record)
    all_sols = sorted(all_sols)
    dedup = [all_sols[i] for i in range(len(all_sols)) if i == 0 or all_sols[i] != all_sols[i-1]]
    return dedup
    #  reps = []
    #  for record in all_sols:
    #      representation = ['.' for w in seq]
    #      for (i, j) in record:
    #          representation[i] = '('
    #          representation[j] = ')'
    #      representation = ''.join(representation)
    #      reps.append(representation)
    #  reps = list(set(reps))
    #  return reps

def to_dot_bracket(seq, record):
    representation = ['.' for w in seq]
    for (i, j) in record:
        representation[i] = '('
        representation[j] = ')'
    representation = ''.join(representation)
    return representation

def visualize_reconstruction(seq, table, record, ell=0):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    pairs = set(record)
    left_p = set([r[0] for r in record])
    right_p = set([r[1] for r in record])
    n = len(seq)
    visited = np.zeros((n, n))
    arrows = [] # quadruples of (x, y, dx, dy)
    v_mark = 1

    stack = [(0,n-1)]
    visited[0, n-1] = v_mark
    while len(stack) > 0:
        v_mark += 1
        i,j = stack.pop()
        #  print(i, j) # color this
        if i + ell >= j:
            continue
        elif (i, j) in pairs:
            # arrow to (i+1, j-1)
            arrows.append((i, j, 1, -1))
            visited[i+1, j-1] = v_mark
            stack.append((i+1, j-1))
        elif i in left_p and j in right_p:
            for k in range(i+1+ell, j):
                if (i, k) in pairs:
                    # draw arrow to (i, k)
                    arrows.append((i, j, 0, k-j))
                    visited[i, k] = v_mark
                    stack.append((i, k))
                    # arrow to (k+1, j)
                    arrows.append((i, j, k+1-i, 0))
                    visited[k+1, j] = v_mark
                    stack.append((k+1, j))
                    break
        elif i in left_p:
            # arrow to (i, j-1)
            arrows.append((i, j, 0, -1))
            visited[i, j-1] = v_mark
            stack.append((i, j-1))
        else:
            # arrow to (i+1, j)
            arrows.append((i, j, 1, 0))
            visited[i+1, j] = v_mark
            stack.append((i+1, j))

    rep = to_dot_bracket(seq, record)
    ax = plt.gca()
    # revert the color for visited entries
    visited[visited > 0] = visited.max() + 1 - visited[visited>0]
    im = ax.imshow(visited, cmap='YlGn')
    threshold = im.norm(im.get_array().max())/2.
    table = np.array(table)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([rep[i] + '\n' + seq[i] + str(i+1) for i in range(n)])
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels([seq[i] + str(i+1) for i in range(n)])
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.set_xticks(np.arange(n+1)-.5, minor=True)
    ax.set_yticks(np.arange(n+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3, zorder=0)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(n):
        for j in range(n):
            text = ax.text(j, i, table[i, j],
                           ha='center', va='center',
                           color='w' if im.norm(visited[i, j])>threshold else 'black')

    for (x, y, dx, dy) in arrows:
        if dy > 0:
            dy -= 0.5
        elif dy < 0:
            dy += 0.5
        if dx > 0:
            dx -= 0.5
        elif dx < 0:
            dx += 0.5
        ax.arrow(y, x, dy, dx, width=0.01, head_width=0.2, alpha=0.5,
                 color='r', zorder=2)
    #  ax.arrow(1, 0, 3, 1, zorder=3)

    plt.show()
    pdb.set_trace()


if __name__ == "__main__":
    GAMMA = set([('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G')])
    string = 'GCAGCAUUCG'
    #  string = 'AUAU'
    #  string = 'GGUCCAC'
    ell = 1
    ans, table = nussinov_table(string, GAMMA, ell)
    print(ans)
    print(table)

    ## single solution
    record = reconstruction(string, table, ell)
    pdb.set_trace()
    print(to_dot_bracket(string, record))

    ## enumerate solutions
    all_solns = reconstruction_all(string, table, GAMMA, ell=ell)
    print('Found {} solutions'.format(len(all_solns)))
    for record in all_solns:
        print(to_dot_bracket(string, record))
        visualize_reconstruction(string, table, record, ell=ell)

