import numpy as np


class DisjointSet:

    def __init__(self, num_vertices):
        self.sizes = np.ones(num_vertices, dtype=int)
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
            return False

        sizes = self.sizes
        parents = self.parents
        if sizes[a] < sizes[b]:
            parents[a] = b
            sizes[b] += sizes[a]
        else:
            parents[b] = a
            sizes[a] += sizes[b]
        return True

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
        unique = np.unique(x)

        # find first occurrences of each element
        indices = {e: len(x) for e in unique}
        for i, e in enumerate(x):
            indices[e] = min(indices[e], i)

        # order elements by frequency, using first occurrences as a tie-breaker
        counts = np.bincount(x)
        ordered = sorted(unique, key=lambda x: (-counts[x], indices[x]))
        assert sorted(ordered) == list(np.unique(x))

        ids = dict([(e, i) for i, e in enumerate(ordered)])
        return np.array([ids[e] for e in x])

    def get_roots(self):
        self._compress()
        return np.unique(self.parents)

    def get_num_components(self):
        return len(self.get_roots())
