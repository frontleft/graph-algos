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
        return len(self.entry_finder) > 0
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
                print(self.cost, G.V, '*')
                self._find_neg_cycle()

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
        print(G)
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


G = EdgeWeightedDigraph(8)
G.add_edge(DirectedEdge(4, 5, 35))
G.add_edge(DirectedEdge(5, 4, -66))
G.add_edge(DirectedEdge(4, 7, 37))
G.add_edge(DirectedEdge(5, 7, 28))
G.add_edge(DirectedEdge(7, 5, 28))
G.add_edge(DirectedEdge(5, 1, 32))
G.add_edge(DirectedEdge(0, 4, 38))
G.add_edge(DirectedEdge(0, 2, 26))
G.add_edge(DirectedEdge(7, 3, 39))
G.add_edge(DirectedEdge(1, 3, 29))
G.add_edge(DirectedEdge(2, 7, 34))
G.add_edge(DirectedEdge(6, 2, 40))
G.add_edge(DirectedEdge(3, 6, 52))
G.add_edge(DirectedEdge(6, 0, 58))
G.add_edge(DirectedEdge(6, 4, 93))


spt = BellmanFordSP(G, 0)
print(spt.has_neg_cycle())
print(spt.get_neg_cycle())
