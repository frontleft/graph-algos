from collections import deque
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
        self._marked = [False for v in range(G.V)]
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


class DirectedCycle():
    def __init__(self, G):
        self.on_stack = [False for v in range(G.V)]
        self.edge_to = [v for v in range(G.V)]
        self.marked = [False for v in range(G.V)]
        self._cycle = []
        for v in range(G.V):
            if not self.marked[v]:
                self._dfs(G, v)

    def _dfs(self, G, v):
        self.on_stack[v] = True
        self.marked[v] = True
        for w in G.adj[v]:
            if self.has_cycle():
                return
            elif not self.marked[w]:
                self.edge_to[w] = v
                self._dfs(G, w)
            elif self.on_stack[w]:
                x = v
                while x != w:
                    self._cycle.append(x)
                    x = self.edge_to[x]
                self._cycle.append(w)
                self._cycle.append(v)
        self.on_stack[v] = False

    def has_cycle(self):
        return len(self._cycle) > 0

    @property
    def cycle(self):
        return self._cycle


class DepthFirstOrder():
    def __init__(self, G):
        self._pre = deque()
        self._post = deque()
        self._reverse_post = deque()
        self.marked = [False for v in range(G.V)]

        for v in range(G.V):
            if not self.marked[v]:
                self._dfs(G, v)

    def _dfs(self, G, v):
        self._pre.appendleft(v)
        self.marked[v] = True
        for w in G.adj[v]:
            if not self.marked[w]:
                self._dfs(G, w)
        self._post.appendleft(v)
        self.reverse_post.append(v)

    @property
    def pre(self):
        return list(self._pre)

    @property
    def post(self):
        return list(self._post)

    @property
    def reverse_post(self):
        return list(self._reverse_post)


class Topological():
    def __init__(self, G):
        self._order = []

        cycle_finder = DirectedCycle(G)
        if not cycle_finder.has_cycle():
            dfs = DepthFirstOrder(G)
            self._order = dfs.reverse_post

    def order(self):
        return self._order

    def is_DAG(self):
        return len(self._order) > 0


class KosarajuSCC():
    ''' Kosaraju's algo for computing strongly connected components '''

    def __init__(self, G):
        self._marked = [False for v in range(G.V)]
        self._id = [v for v in range(G.V)]
        self._count = 0
        order = DepthFirstOrder(G.reverse())
        for s in order.reverse_post:
            if not self._marked[s]:
                self._dfs(G, s)
                self._count += 1

    def _dfs(self, G, v):
        self._marked[v] = True
        self._id[v] = self._count
        for w in G.adj[v]:
            if not self._marked[w]:
                self._dfs(G, w)

    def strongly_connected(self, v, w):
        return self._id[v] == self._id[w]

    def identifier(self, v):
        return self._id[v]

    @property
    def count(self):
        return self._count


class Hamiltonian():
    ''' Tests for Hamiltonian Path by taking the topological order and checking
    to see if edges exist for each consecutive pair 
    '''

    def __init__(self, G):
        self.path = True
        topo = Topological(G)
        if topo.is_DAG():
            order = topo.order()
            for idx, v in enumerate(order):
                if (idx == len(order) - 1):
                    continue
                if order[idx+1] not in G.adj(v):
                    self.path = False
        else:
            self.path = False

    def exist(self):
        return self.path
