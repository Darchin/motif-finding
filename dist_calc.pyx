cpdef int distance(str s1, str s2):
    cdef int dist = 0
    cdef int n = len(s1)
    for i in range(n):
        if s1[i] != s2[i]: dist += 1
    return dist