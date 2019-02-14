import numpy as np


class DisjointSet:

    def __init__(self, num_vertices):

        self.ranks = np.zeros(num_vertices).astype(np.int)
        self.parents = np.arange(num_vertices)

    def find(self, index):

        parents = self.parents
        parent = parents[index]
        while parent != parents[parent]:
            parent = parents[parent]
        parents[index] = parent
        return parent

    def merge(self, a, b):

        a = self.find(a)
        b = self.find(b)
        if a == b:
            return

        ranks = self.ranks
        parents = self.parents

        if ranks[a] < ranks[b]:
            parents[a] = b
        elif ranks[a] > ranks[b]:
            parents[b] = a
        else:
            parents[b] = a
            ranks[a] += 1

    def _compress(self):

        a = self.parents
        b = a[a]
        while (a != b).any():
            a = b
            b = a[a]
        self.parents = a

    def get_components(self, relabel=False):

        self._compress()
        if not relabel:
            return self.parents

        x = np.copy(self.parents)
        ids = dict([(e, i) for i, e in enumerate(np.unique(x))])
        for i in range(len(x)):
            x[i] = ids[x[i]]
        return x

    def get_roots(self):

        self._compress()
        return np.unique(self.parents)

    def get_num_components(self):

        return len(self.get_roots())
