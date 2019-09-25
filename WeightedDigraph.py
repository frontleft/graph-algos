import heapq
import itertools
from collections import deque
''' Python3 implementation of weighted graph algos'''


class DirectedEdge():
    def __init__(self, v, w, weight):
        self.v = v
        self.w = w
        self.weight = weight

    def __str__(self):
        return f'{self.v}->{self.w}, {self.weight}'

    def other(self, vertex):
        if (vertex == self.v):
            return self.w
        elif (vertex == self.w):
            return self.v
        else:
            raise IndexError()

    def either(self):
        return self.v

    @property
    def from_v(self):
        return self.v

    @property
    def to_v(self):
        return self.w


class EdgeWeightedDigraph():
    def __init__(self, V):
        self._V = V  # Number of vertices in the graph
        self._E = 0  # Number of edges in the graph
        self.adj = [[] for _ in range(0, V)]  # Adjacency list setup

    def __str__(self):
        s = f'{self.V} vertices, {self.E} edges \n'
        for v in range(0, self.V):
            s += f'{v}: '
            for edge in self.adj[v]:
                s += f' ({edge.to_v}, w{edge.weight}) '
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

    def add_edge(self, edge):
        '''
        Adds a weighted directed edge v->w to the digraph

        Parameters:
        edge    Edge    edge to be added

        Return:
        void
        '''
        self.adj[edge.from_v].append(edge)
        self.E += 1

    def get_edges(self):
        edges = []
        for v in range(self.V):
            for edge in self.adj[v]:
                edges.append(edge)
        return edges


class MinPriorityQueue():
    REMOVED = '<removed-task>'

    def __init__(self):
        self.heap = []
        self.entry_finder = {}
        self.counter = itertools.count()

    def remove_task(self, task):
        '''Mark existing task as removed; raise keyError if not found'''
        entry = self.entry_finder.pop(task)
        entry[-1] = self.REMOVED

    def insert_task(self, task, priority=0):
        ''' add a new task or update priority of an existing task'''
        if task in self.entry_finder:
            self.remove_task(task)
        count = next(self.counter)
        entry = [priority, count, task]
        self.entry_finder[task] = entry
        heapq.heappush(self.heap, entry)

    def del_min(self):
        ''' remove & return lowest priority task; keyError if empty'''
        while self.heap:
            priority, count, task = heapq.heappop(self.heap)
            if task is not self.REMOVED:
                del self.entry_finder[task]
                return task
        raise KeyError('pop from empty PQ')

    def is_empty(self):
        ''' determines if there are tasks in the PQ based on entry finder'''
        return len(self.entry_finder) <= 0
        # while self.heap[0][-1] == self.REMOVED:
        #     heapq.heappop(self.heap)
        # if len(self.heap) == 0:
        #     return True
        # return False

    def contians(self, task):
        if task in self.entry_finder:
            return True
        return False


class DijkstraSP():
    ''' Dijkstra's SP algo implemented w/ priority queue for O(ELogV) complexity '''

    def __init__(self, G, s):
        self.edge_to = [None for v in range(G.V)]
        self.dist_to = [float('inf') for v in range(G.V)]
        self.dist_to[s] = 0
        self.pq = MinPriorityQueue()

        self.pq.insert_task(s, 0)
        while not self.pq.is_empty():

            self._relax(G, self.pq.del_min())

    def _relax(self, G, v):
        for edge in G.adj[v]:
            w = edge.to_v
            if self.dist_to[w] > self.dist_to[v] + edge.weight:
                self.dist_to[w] = self.dist_to[v] + edge.weight
                self.edge_to[w] = edge
                self.pq.insert_task(w, self.dist_to[w])

    def has_path_to(self, v):
        return self.dist_to[v] < float('inf')

    def path_to(self, v):
        if not self.has_path_to(v):
            return None
        path = []
        e = self.edge_to[v]
        while e != None:
            path.append(e)
            e = self.edge_to[e.from_v]
        return path


class BellmanFordSP():
    '''Bellman-Ford algorithm w/ a queue to find shortest path in a edge-weighted digraph with negative edges
    O(EV) complexity
    '''

    def __init__(self, G, s):
        self.edge_to = [None for v in range(G.V)]
        self.dist_to = [float('inf') for v in range(G.V)]
        self.dist_to[s] = 0
        self.on_q = [False for v in range(G.V)]
        self.q = deque()
        self.cost = 0  # num of times relax has been called
        self.cycle = []
        # initialize the state of the queue
        self.q.append(s)
        self.on_q[s] = True

        while len(self.q) > 0 and not self.has_neg_cycle():
            v = self.q.popleft()

            self.on_q[v] = False
            self.relax(G, v)

    def has_neg_cycle(self):
        return len(self.cycle) > 0

    def get_neg_cycle(self):
        return self.cycle

    def _find_neg_cycle(self):
        ''' builds a new graph of just the edges on the shortest paths, searches for cycle
        if there is a cycle, it must be negative weight'''
        V = len(self.edge_to)
        spt = EdgeWeightedDigraph(V)
        for v in range(V):
            if self.edge_to[v]:
                spt.add_edge(self.edge_to[v])
        cf = EdgeWeightedCycleFinder(spt)
        self.cycle = cf.cycle

    def relax(self, G, v):
        self.cost += 1
        for edge in G.adj[v]:
            w = edge.to_v
            if self.dist_to[w] > self.dist_to[v] + edge.weight:
                self.dist_to[w] = self.dist_to[v] + edge.weight
                self.edge_to[w] = edge
                if not self.on_q[w]:
                    self.q.append(w)
                    self.on_q[w] = True
            if self.cost % G.V == 0:
                self._find_neg_cycle()

    def has_path_to(self, v):
        return self.dist_to[v] < float('inf')

    def path_to(self, v):
        if not self.has_path_to(v) or self.has_neg_cycle():
            return None
        path = []
        e = self.edge_to[v]
        while e != None:
            path.append(e)
            e = self.edge_to[e.from_v]
        return path


class EdgeWeightedCycleFinder():
    def __init__(self, G):
        self.on_stack = [False for v in range(G.V)]
        self.edge_to = [None for v in range(G.V)]
        self.marked = [False for v in range(G.V)]
        self._cycle = []
        for v in range(G.V):
            if not self.marked[v]:
                self._dfs(G, v)

    def _dfs(self, G, v):
        self.on_stack[v] = True
        self.marked[v] = True
        for edge in G.adj[v]:
            w = edge.to_v
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


class LazyPrimsMST():
    '''
    ** INVALID FOR Directed Graphs; Only valid for undirected graphs **
    '''

    def __init__(self, G):
        self._marked = [False for v in range(G.V)]
        self.mst = deque()
        self.q = MinPriorityQueue()

        # initialize the queue (arbitrarily w/ vertex 0)
        self.visit(G, 0)
        while not self.q.is_empty():
            e = self.q.del_min()
            v = e.either()
            w = e.other(v)
            if not self._marked[v] or not self._marked[w]:
                self.mst.append(e)  # add to the MST
            if not self._marked[v]:
                self.visit(G, v)
            if not self._marked[w]:
                self.visit(G, w)

    @property
    def weight(self):
        '''
        Returns the weight of the edges in the minimum spanning tree
        '''
        w = 0
        for edge in list(self.mst):
            w += edge.weight
        return w

    @property
    def edges(self):
        return self.mst

    def visit(self, G, v):
        '''
        Visits a vertex, adds it to the MST, and enqueues its valid  adjacent edges on the PQ

        Parameters:
        G       EdgeWeightedDigraph         Edge-weighted digraph from which to construct the MST.
        v       Int                         Index of the vertex to visit.

        Return
        void
        '''
        self._marked[v] = True
        for edge in G.adj[v]:
            if not self._marked[edge.other(v)]:
                self.q.insert_task(edge, edge.weight)


class EagerPrimsMST():
    '''
    ** INVALID FOR Directed Graphs; Only valid for undirected graphs **
    '''

    def __init__(self, G):
        self.marked = [False for v in range(G.V)]
        self.edge_to = [None for v in range(G.V)]
        self.dist_to = [float('inf') for v in range(G.V)]
        self.q = MinPriorityQueue()
        self._weight = 0

        # initialize the PQ (arbitrarily w/ vertex 0)
        self.dist_to[0] = 0
        self.q.insert_task(0, 0)

        while not self.q.is_empty():
            self.visit(G, self.q.del_min())

    def visit(self, G, v):
        self.marked[v] = True
        for edge in G.adj[v]:
            w = edge.other()
            if not self.marked[w] and edge.weight < self.dist_to[w]:
                # eager update of weight
                if self.dist_to[w] == float('inf'):
                    self._weight += edge.weight
                else:
                    self._weight -= (self.dist_to[w] - edge.weight)

                self.edge_to[w] = edge
                self.dist_to[w] = edge.weight
                self.q.insert_task(w, self.dist_to[w])

    @property
    def edges(self):
        return self.edge_to[1:]

    @property
    def weight(self):
        return self._weight


class UnionFind():
    def __init__(self, N):
        self.id = [i for i in range(N)]
        self.weight = [1 for i in range(N)]
        self._count = N

    @property
    def count(self):
        return self._count

    def union(self, u, v):
        u_parent = self.find(u)
        v_parent = self.find(v)
        if self.weight[u_parent] > self.weight[v_parent]:
            self.id[v_parent] = u_parent
            self.weight[u_parent] += self.weight[v_parent]
        else:
            self.id[u_parent] = v_parent
            self.weight[v_parent] += self.weight[u_parent]
        self.count -= 1

    def find(self, u):
        while u != self.id[u]:
            u = self.id[u]
        return u

    def connected(self, u, v):
        # if they have the same parent, true; else false
        return self.find(u) == self.find(v)


class KruskalsMST():
    '''
    ** INVALID FOR Directed Graphs; Only valid for undirected graphs **
    '''

    def __init__(self, G):
        self.mst = []
        self.q = MinPriorityQueue()
        self.uf = UnionFind(G.V)
        self._weight = 0

        while not self.q.is_empty():
            edge = self.q.del_min()
            v = edge.either()
            w = edge.other(v)
            if self.uf.connected(v, w):
                continue
            self.uf.union(v, w)
            self.mst.append(edge)
            self._weight += edge.weight

    @property
    def edges(self):
        return self.mst

    @property
    def weight(self):
        return self._weight


class FloydWarshall():
    def __init__(self, G):
        self.dist = [[float('inf') for i in range(G.V)] for j in range(G.V)]
        self.next = [[None for i in range(G.V)] for j in range(G.V)]

        # * Initialize all distances to their single-edge weights in the adj. matrix
        for edge in G.get_edges():
            self.dist[edge.from_v][edge.to_v] = edge.weight
            self.next[edge.from_v][edge.to_v] = edge.to_v
        # * Initialize all self-loops to 0
        for v in range(G.V):
            self.dist[v][v] = 0
            self.next[v][v] = v

        for k in range(1, G.V):
            for i in range(1, G.V):
                for j in range(1, G.V):
                    if self.dist[i][j] > self.dist[i][k] + self.dist[k][j]:
                        self.dist[i][j] = self.dist[i][k] + self.dist[k][j]
                        self.next[i][j] = self.next[i][k]

    def distance(self, u, v):
        return self.dist[u][v]

    def path(self, u, v):
        if self.next[u][v] == None:
            return []
        output = [u]
        while u != v:
            u = self.next[u][v]
            output.append(u)
        return output
