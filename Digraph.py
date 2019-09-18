'''
Python3 implementation of unweighted directed graph data type
'''


class Digraph():
    def __init__(self, V):
        self._V = V  # Number of vertices in the graph
        self._E = 0  # Number of edges in the graph
        self.adj = [[] for _ in range(0, V)]  # Adjacency list setup

    def __str__(self):
        s = f'{self.V} vertices, {self.E} edges \n'
        for v in range(0, self.V):
            s += f'{v}: '
            for w in self.adj[v]:
                s += f'{w} '
            s += '\n'
        return s

    @property
    def V(self):
        return self._V

    @V.setter
    def V(self, V):
        self._V = V

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, E):
        self._E = E

    def add_edge(self, v, w):
        '''
        Adds a directed edge v->w to the digraph

        Parameters:
        v       Int     integer representing from vertex
        w       Int     integer representing to vertex

        Return:
        void
        '''
        self.adj[v].append(w)
        self.E += 1

    def reverse(self):
        '''
        Returns a copy of the graph with all edges reversed

        Return:
        Digraph            Reversed digraph
        '''
        R = Digraph(self.V)
        for v in range(0, self.V):
            for w in self.adj[v]:
                R.add_edge(w, v)
        return R


class DirectedDFS():
    def __init__(self, G, sources):
        self.G = G  # A Digraph
        self.sources = sources  # A list of sources
        self._marked = [False for _ in range(G.V)]
        for s in sources:
            if not self.is_marked(s):
                self._dfs(G, s)
        for v in range(self.G.V):
            if self.is_marked(v):
                print(f'{v} ')

    def _dfs(self, G, v):
        self._marked[v] = True
        for w in G.adj[v]:
            if not self._marked[w]:
                self._dfs(G, w)

    def is_marked(self, v):
        return self._marked[v]
