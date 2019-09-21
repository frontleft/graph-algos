from collections import deque
'''
Python3 Implementation of flow network algorithms (Ford Fulkerson, etc.)
'''


class FlowEdge():
    def __init__(self, v, w, capacity, flow):
        self._v = v
        self._w = w
        self._capacity = capacity
        self._flow = flow

    def __str__(self):
        return f'{self._v}->{self._w} ({self._flow}/{self._capacity})'

    @property
    def v(self):
        return self._v

    @property
    def from_v(self):
        return self._v

    @property
    def w(self):
        return self._w

    def to_v(self):
        return self._w

    @property
    def capacity(self):
        return self._capacity

    @capacity.setter
    def capacity(self, new_capacity):
        if (new_capacity >= 0):
            self._capacity = new_capacity
        else:
            raise ValueError('capacity must be >= 0')

    @property
    def flow(self):
        return self._flow

    @flow.setter
    def flow(self, new_flow):
        if (new_flow >= 0):
            self._flow = new_flow
        else:
            raise ValueError('flow must be >= 0')

    def other(self, vertex):
        if (vertex == self._v):
            return self._w
        elif (vertex == self._w):
            return self._v
        else:
            raise IndexError('Invalid endpoint')

    def residual_capacity_to(self, vertex):
        if (vertex == self._v):
            return self._flow
        elif (vertex == self._w):
            return self._capacity - self._flow
        else:
            raise IndexError('Invalid endpoint')

    def add_resid_flow_to(self, vertex, delta):
        if (vertex == self._v):
            self.flow -= delta
        elif (vertex == self._w):
            self.flow += delta
        else:
            raise IndexError('Invalid endpoint')


class FlowNetwork():
    def __init__(self, V):
        self._V = V
        self._E = 0
        self._adj = [[] for v in range(V)]

    def __str__(self):
        s = f''
        return s

    @property
    def V(self):
        return self._V

    @property
    def E(self):
        return self._V

    @E.setter
    def E(self, new_E):
        self._E = new_E

    @property
    def adj(self):
        return self._adj

    def add_edge(self, edge):
        self.E += 1
        v = edge.from_v
        w = edge.to_v
        self.adj[v].appened(edge)
        self.adj[w].append(edge)


class FordFulkerson():
    def __init__(self, G, s, t):
        self._value = 0
        self._marked = [False for v in range(G.V)]
        self.edge_to = [None for v in range(G.V)]

        while self.has_augmenting_path(G, s, t):
            bottle = float('inf')
            v = t
            while v != s:  # compute bottleneck capacity
                bottle = min(bottle, self.edge_to[v].residual_capacity_to(v))
                v = self.edge_to[v].other(v)

            v = t
            while v != s:  # augment the flow
                self.edge_to[v].add_resid_flow_to(v, bottle)
                v = self.edge_to[v].other(v)

            self._value += bottle

    def has_augmenting_path(self, G, s, t):
        edge_to = [None for v in range(G.V)]
        marked = [False for v in range(G.V)]

        q = deque()
        q.append(s)
        marked[s] = True
        while len(q) > 0:
            v = q.popleft()
            for edge in G.adj[v]:
                w = edge.other(v)
                # while there is a path from s->w
                if edge.residual_capacity_to(w) > 0 and not marked[w]:
                    edge_to[w] = edge
                    marked[w] = True
                    q.append(w)
        return marked[t]

    def in_cut(self, v):
        return self._marked[v]

    @property
    def value(self):
        return self._value
