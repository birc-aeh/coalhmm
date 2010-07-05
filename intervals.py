def do_on_all_distributions(S, I, visit):
    buffer = [S[0]] * I
    def f(S, S_i, I, I_i):
        if S_i == len(S):
            return
        for i in xrange(I_i, I):
            r = buffer[i:]
            buffer[i:] = [S[S_i]] * len(r)
            visit(buffer)
            f(S, S_i+1, I, i+1)
            buffer[i:] = r
    visit(buffer)
    for i in xrange(1,len(S)):
        f(S, i, I, 0)

if __name__ == "__main__":
    def p(s):
        print s
    do_on_all_distributions("abc", 4, p)
