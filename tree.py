iset = frozenset

def make_tree(G, s, side):
    """ Given a path, s, through the graph, G, create a tree with timestamps.
        The tree is either a leaf (an iset with 1 member) or a node: a tuple
        with the timestamp and an iset of children.
    """
    tree = None
    initial = G.project_state(0, G.initial(0), side)
    used = set()
    for i in xrange(len(s)):
        ce, cs = s[i]
        pe, ps = s[i-1]
        if i == 0:
            B = G.projected(ce, cs, side)
            if len(initial) != len(B):
                joined_from_the_start = initial - B
                tree = (0, joined_from_the_start)
                for x in joined_from_the_start:
                    used.add(x)
        elif cs != ps:
            A = G.projected(pe, ps, side)
            B = G.projected(ce, cs, side)
            A = iset([x for x in A if len(x) == 1 and x not in used])
            joined = A-B
            if len(joined) <= 1:
                continue
            for x in joined:
                used.add(x)
            if tree == None:
                tree = (i, joined)
            else:
                tree = (i, iset([tree]).union(joined))
    if len(initial) == len(used):
        return tree
    rest_joined_at = len(s)
    rest = iset([x for x in initial if not x in used])
    return (rest_joined_at, (tree and iset([tree]) or iset()).union(rest))

def tree_to_newick(t):
    """Converts a tree into newick format """
    if isinstance(t, iset):
        return "".join([str(x) for x in t])
    else:
        return "(" + ", ".join(map(tree_to_newick, t[1])) + "):" + str(t[0])

if __name__ == "__main__":
    from model import *
    theta = 2*30000.0 * 25 * 1e-9
    C = 1.0 / theta
    R = 1.5e-8 / 1.0e-9

    _, G = build_epoch_seperated_scc(3, [[0, 1, 1], [0, 0, 0]])
    path = ((0, 0), (1, 2), (1, 2), (1, 2), (2, 24), (2, 15), (2, 15))
    assert "((0, 1):5, 2):7" == tree_to_newick(make_tree(G, path, 0))
    assert "(0, (1, 2):1):5" == tree_to_newick(make_tree(G, path, 1))

