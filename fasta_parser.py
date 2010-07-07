def readAlignment(filename):
    '''Read a multiple sequence alignment from a fasta file, storing it as
    a dictionary that maps the species name to the sequence.'''

    results = dict()
    name = None
    buffer = []
    for x in (open(filename)).readlines():
        if x.startswith(">"):
            if name != None:
                results[name] = "".join(buffer)
            name = x[1:].strip()
            buffer = []
        else:
            buffer.append(x.strip())
    results[name] = "".join(buffer)

    return results


