from collections import deque
'''
Python3 Implementation of flow network algorithms (Ford Fulkerson, etc.)
'''


class FlowEdge():
    def __init__(self, v, w, capacity, flow=0):
        '''
        Parameters:
        v           Int     Index of 'from' vertex
        w           Int     Index of 'to' vertex
        capacity    Int     Forward capacity of the edge
        flow        Int     Amount of forward flow currently going through the edge
        '''
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

    @property
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
        '''
        Returns the index of the other vertex on the edge (i.e. if edge v->w given v returns w)

        Parameters:
        vertex      Int     Index of a vertex on the edge

        Return:
        Int                 Index of the other vertex
        '''
        if (vertex == self._v):
            return self._w
        elif (vertex == self._w):
            return self._v
        else:
            raise IndexError('Invalid endpoint')

    def residual_capacity_to(self, vertex):
        '''
        Returns the amount of residual capacity (i.e. capacity - flow) to the provided vertex

        Parameters:
        vertex      Int     Index of a vertex on the edge

        Return:
        Int                 Amount of flow that can be directed to the vertex
        '''
        if (vertex == self._v):
            return self._flow
        elif (vertex == self._w):
            return self._capacity - self._flow
        else:
            raise IndexError('Invalid endpoint')

    def add_resid_flow_to(self, vertex, delta):
        '''
        Adds the 'delta' amount of flow along the edge in the direction of the vertex

        Parameters:
        vertex      Int     Index of a vertex on the edge
        delta       Int     if > 0, Amount of flow to direct to the vertex, if < 0, amount of flow to remove from the vertex

        Return:
        void
        '''
        if (vertex == self._v):
            self.flow -= delta
        elif (vertex == self._w):
            self.flow += delta
        else:
            raise IndexError('Invalid endpoint')


class FlowNetwork():
    def __init__(self, V):
        '''
        Parameters:
        V      Int      Number of vertices within the flow network
        '''
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
        '''
        Add an edge to the flow network (adds to adj. list of both vertices)

        Parameters:
        edge        FlowEdge        Flow edge to be added to the network

        Return:
        void
        '''
        self.E += 1
        v = edge.from_v
        w = edge.to_v

        self.adj[v].append(edge)
        self.adj[w].append(edge)


class FordFulkerson():
    def __init__(self, G, s, t):
        '''
        Parameters:
        G       FlowNetwork     Flow network to search for augmenting path within
        s       Int             Source vertex
        t       Int             Sink vertex  
        '''
        self._value = 0
        self._marked = [False for v in range(G.V)]
        self.edge_to = [None for v in range(G.V)]
        self._G = G
        self._s = s
        self._t = t

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

    def __str__(self):
        s = f'Max flow from {self._s} to {self._t}: \n'
        for v in range(self._G.V):
            for edge in self._G.adj[v]:
                if v == edge.from_v and edge.flow > 0:
                    s += str(edge) + '\n'
        s += f'Max flow value = {self._value}\n'
        return s

    def has_augmenting_path(self, G, s, t):
        '''
        Parameters:
        G       FlowNetwork     Flow network to search for augmenting path within
        s       Int             Source vertex
        t       Int             Sink vertex  

        Return:
        bool                    True if there is a augmenting path from s->t else false
        '''
        self.edge_to = [None for v in range(G.V)]
        self._marked = [False for v in range(G.V)]

        q = deque()
        q.append(s)
        self._marked[s] = True
        while len(q) > 0:
            v = q.popleft()
            for edge in G.adj[v]:
                w = edge.other(v)
                # while there is a path from s->w, save the edge, mark it as known & add to the q
                if edge.residual_capacity_to(w) > 0 and self._marked[w] == False:
                    self.edge_to[w] = edge
                    self._marked[w] = True
                    q.append(w)
        return self._marked[t]

    def in_cut(self, v):
        return self._marked[v]

    @property
    def value(self):
        return self._value


G = FlowNetwork(6)
G.add_edge(FlowEdge(0, 1, 2))
G.add_edge(FlowEdge(0, 2, 3))
G.add_edge(FlowEdge(1, 3, 3))
G.add_edge(FlowEdge(1, 4, 1))
G.add_edge(FlowEdge(2, 3, 1))
G.add_edge(FlowEdge(2, 4, 1))
G.add_edge(FlowEdge(3, 5, 2))
G.add_edge(FlowEdge(4, 5, 3))
ff_g = FordFulkerson(G, 0, 5)
print(ff_g)
