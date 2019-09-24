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

        # mutates the edge_to array to reflect an augmenting path
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
        Mutates the edge_to instance variable to reflect an augmenting path if it exists. This method uses
        BFS to find augmenting paths (i.e. Edmonds-Karp O(E^2 * V)).

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

    def min_cut(self):
        '''
        Finds the minimum cut in a flow network by first finding the max flow. When FF terminates it has found a max flow and the last BFS will have marked the reachable vertices in the residual graph
        '''
        cut = []
        for v, marked in enumerate(self._marked):
            if marked:
                cut.append(v)
        return cut

    def in_cut(self, v):
        return self._marked[v]

    @property
    def value(self):
        return self._value


class Dinics():
    ''' Constructs layered graph from residual graph via BFS. Compute blocking flow in GF, augment f by g. Blocking flow means there is no path in G such that there is room left on every path. O(V^2 * E)
    '''

    def __init__(self, G, s, t):
        self.dist_to = [float('inf') if v != 0 else 0 for v in range(G.V)]
        self.total_flow = 0
        self._s = s
        self._t = t
        self._G = G

        # assumes all edges' flow set to 0

        #   construct level graph (Gl)from residual graph of G
        while self.construct_level(G, s, t):
            #   find a blocking flow (f') in Gl & augment flow
            flow = self.augment_flow(G, s, t, float('inf'))
            while (flow > 0):
                self.total_flow += flow
                flow = self.augment_flow(G, s, t, float('inf'))

    def __str__(self):
        s = f'Max flow from {self._s} to {self._t}: \n'
        for v in range(self._G.V):
            for edge in self._G.adj[v]:
                if v == edge.from_v and edge.flow > 0:
                    s += str(edge) + '\n'
        s += f'Max flow value = {self.total_flow}\n'
        return s

    def construct_level(self, G, s, t):
        '''
        Calculates distances in the level graph; returns True if there is an s-t path in the graph, false otherwise
        '''
        self.dist_to = [float('inf') if v != 0 else 0 for v in range(G.V)]
        q = deque()  # initialize a queue for the BFS
        q.append(s)
        while len(q) > 0:
            v = q.popleft()
            for edge in G.adj[v]:
                w = edge.other(v)
                if self.dist_to[w] == float('inf') and edge.residual_capacity_to(w) > 0:
                    self.dist_to[w] = self.dist_to[v] + 1
                    q.append(w)
        return self.dist_to[t] != float('inf')

    def augment_flow(self, G, u, t, flow):
        '''
        Conducts a DFS through the residual graph
        '''
        # base case
        if u == t:
            return flow

        for edge in G.adj[u]:
            w = edge.other(u)
            residual_capacity = edge.residual_capacity_to(w)
            if self.dist_to[w] == self.dist_to[u] + 1 and residual_capacity > 0:
                cur_flow = min(flow, residual_capacity)
                temp_flow = self.augment_flow(G, w, t, cur_flow)

                if temp_flow > 0:
                    edge.add_resid_flow_to(w, temp_flow)
                    return temp_flow

        return 0


class DoublyLinkedList():
    def __init__(self):
        self.head = None
        self.tail = None

    def __str__(self):
        s = ''
        cur = self.head
        while cur != None:
            s += f'- {cur}'
            cur = cur.next
        return s

    def is_empty(self):
        return self.head == None

    def insert_tail(self, item):
        new_node = DoublyLinkedListNode(item)
        if not self.tail:
            self.head = new_node
            self.tail = new_node
        else:
            new_node.prev = self.tail
            self.tail.next = new_node
            self.tail = new_node
        return new_node

    def insert_head(self, item):
        new_node = DoublyLinkedListNode(item)
        if not self.head:
            self.head = new_node
            self.tail = new_node
        else:
            new_node.next = self.head
            self.head.prev = new_node
            self.head = new_node
        return new_node

    def peek_head(self):
        return self.head

    def pop_head(self):
        if not self.head:
            raise IndexError('Cannot pop from empty list')
        tmp = self.head
        if tmp.next:
            self.head = tmp.next
            self.head.prev = None
            tmp.next = None
        else:
            self.tail = None
            self.head = None
        return tmp

    def pop_tail(self):
        if not self.tail:
            raise IndexError('Cannot pop from empty list')
        tmp = self.tail
        if tmp.prev:
            self.tail = tmp.prev
            self.tail.next = None
            tmp.prev = None
        else:
            self.tail = None
            self.head = None
        return tmp

    def remove(self, node):
        if node == None:
            return None
        if not self.head:
            raise IndexError('Cannot remove from empty list')
        if node == self.head and node == self.tail:
            self.head = None
            self.tail = None
        elif node == self.head:
            self.head = node.next
        elif node == self.tail:
            self.tail = node.prev
        else:
            next_node = node.next
            prev_node = node.prev
            next_node.prev = prev_node
            prev_node.next = next_node
        node.prev = None
        node.next = None
        return node


class DoublyLinkedListNode():
    def __init__(self, item=None):
        self.next = None
        self.prev = None
        self.item = item

    def __str__(self):
        return f'({self.item})'


class PushRelabel():
    '''
    Computes the maximum flow in a flow network using push-relabel / preflow-push algorithm in O(V^2*E)
    '''

    def __init__(self, G, s, t):
        '''
        Assumes the flow network begins with flow of 0

        Invariant:
        - height(source) = G.V at all times
        - height(sink) = 0
        - each edge (v,w) of the residual network (with + resid capacity) has a height
          differential of <= +1 

        Parameters
        G           FlowNetwork     The flow network to compute the max flow within
        s           Int             The index of the source vertex
        t           Int             The index of the sink vertex
        '''
        self.excess_flows = [0 for v in range(G.V)]
        self.heights = [0 if v != s else G.V for v in range(G.V)]
        self.overflow_list = DoublyLinkedList()
        self.overflows = [None for v in range(G.V)]
        self._G = G
        self._s = s
        self._t = t

        # initialize the edges coming out of the source with a preflow equivalent to their capacity
        for edge in G.adj[s]:
            w = edge.other(s)
            preflow = edge.residual_capacity_to(w)
            edge.add_resid_flow_to(w, preflow)
            self.excess_flows[w] += preflow
            self.excess_flows[s] -= preflow
            # source initially has negative excess flow
            self.overflows[w] = self.overflow_list.insert_tail(w)

        while not self.overflow_list.is_empty():
            # take vertex from front of overflow list (repeatedly; pops when excess is removed in push)
            u = self.overflow_list.peek_head().item
            pushed = False
            # attempt to push if u has an outgoing edge that is downhill & has capacity
            for edge in G.adj[u]:
                w = edge.other(u)
                if edge.residual_capacity_to(w) > 0 and self.heights[u] == self.heights[w] + 1:
                    self.push(u, edge)
                    pushed = True
                    break
            if not pushed:  # else relabel u
                if u == self._s or u == self._t:  # can't change height of the sink or source
                    continue
                self.relabel(G, u)

    def __str__(self):
        s = f'Max flow from {self._s} to {self._t}: \n'
        for v in range(self._G.V):
            for edge in self._G.adj[v]:
                if v == edge.from_v and edge.flow > 0:
                    s += str(edge) + '\n'
        return s

    def push(self, u, edge):
        v = edge.other(u)
        delta_flow = min(self.excess_flows[u], edge.residual_capacity_to(v))
        edge.add_resid_flow_to(v, delta_flow)

        self.excess_flows[u] -= delta_flow
        # * if u is now not overflowing, need to remove from the overflow list
        if self.excess_flows[u] == 0:
            node = self.overflows[u]
            self.overflow_list.remove(node)
            self.overflows[u] = None

        # * if not already in the overflow list, need to add v & save pointer
        self.excess_flows[v] += delta_flow
        if self.excess_flows[v] > 0 and not self.overflows[v] and (v != self._s and v != self._t):
            self.overflows[v] = self.overflow_list.insert_tail(v)

    def relabel(self, G, u):
        min_height = self.heights[u]
        for edge in G.adj[u]:
            v = edge.other(u)
            if edge.residual_capacity_to(v) > 0:
                min_height = min(min_height, self.heights[v])
        self.heights[u] = min_height + 1


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

G2 = FlowNetwork(6)
G2.add_edge(FlowEdge(0, 1, 2))
G2.add_edge(FlowEdge(0, 2, 3))
G2.add_edge(FlowEdge(1, 3, 3))
G2.add_edge(FlowEdge(1, 4, 1))
G2.add_edge(FlowEdge(2, 3, 1))
G2.add_edge(FlowEdge(2, 4, 1))
G2.add_edge(FlowEdge(3, 5, 2))
G2.add_edge(FlowEdge(4, 5, 3))
d_g = Dinics(G2, 0, 5)
print(d_g)

G3 = FlowNetwork(6)
G3.add_edge(FlowEdge(0, 1, 16))
G3.add_edge(FlowEdge(0, 2, 13))
G3.add_edge(FlowEdge(1, 2, 10))
G3.add_edge(FlowEdge(1, 3, 12))
G3.add_edge(FlowEdge(2, 1, 4))
G3.add_edge(FlowEdge(2, 4, 14))
G3.add_edge(FlowEdge(3, 2, 9))
G3.add_edge(FlowEdge(3, 5, 20))
G3.add_edge(FlowEdge(4, 3, 7))
G3.add_edge(FlowEdge(4, 5, 4))
dg2 = Dinics(G3, 0, 5)
print(dg2)

G4 = FlowNetwork(6)
G4.add_edge(FlowEdge(0, 1, 16))
G4.add_edge(FlowEdge(0, 2, 13))
G4.add_edge(FlowEdge(1, 2, 10))
G4.add_edge(FlowEdge(1, 3, 12))
G4.add_edge(FlowEdge(2, 1, 4))
G4.add_edge(FlowEdge(2, 4, 14))
G4.add_edge(FlowEdge(3, 2, 9))
G4.add_edge(FlowEdge(3, 5, 20))
G4.add_edge(FlowEdge(4, 3, 7))
G4.add_edge(FlowEdge(4, 5, 4))
p_p = PushRelabel(G4, 0, 5)
print(p_p)
