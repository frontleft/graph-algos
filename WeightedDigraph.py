import heapq
import itertools
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
