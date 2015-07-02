@profile
def func():
    a = []
    b = []
    for i in xrange(500):
        a.append(range(i))
        b.extend(range(i))
    del b
    return a
