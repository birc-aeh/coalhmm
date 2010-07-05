asd = []
def print_all_distributions(S, I):
    def f(S, S_i, I, I_i):
        print asd
        if S_i == len(S):
            return
        for i in xrange(I_i, I):
            asd.append((i, S[S_i]))
            f(S, S_i+1, I, i+1)
            asd.pop()
    def g(S, S_i, I, I_i):
        print asd
        if S_i == len(S):
            return
        for i in xrange(I_i, I):
            r = asd[i:]
            asd[i:] = [S[S_i]] * len(r)
            g(S, S_i+1, I, i+1)
            asd[i:] = r
    asd = [S[0]] * I
    for i in xrange(I):
        g(S, i, I, 0)

print_all_distributions("abc", 3)
