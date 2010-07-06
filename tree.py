from sets import ImmutableSet as iset

def make_tree(G, s, side):
    """ Given a path, s, through the graph, G, create a tree with timestamps.
        The tree is either a leaf (an iset with 1 member) or a node: a tuple
        with the timestamp and an iset of children.
    """
    tree = None
    initial = G.project_state(G.initial(), side)
    used = set()
    for i in xrange(len(s)):
        if i == 0:
            B = G.projected(s[i], side)
            if len(initial) != len(B):
                joined_from_the_start = initial - B
                tree = (0, joined_from_the_start)
                for x in joined_from_the_start:
                    used.add(x)
        elif s[i] != s[i-1]:
            A = G.projected(s[i-1], side)
            B = G.projected(s[i], side)
            A = iset([x for x in A if len(x) == 1])
            joined = A-B
            if len(joined) == 0:
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
