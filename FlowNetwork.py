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
    def w(self):
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

    def residual_capacity(self, vertex):
        pass

    def add_resid_flow_to(self, vertex, delta):
        pass


class FlowNetwork():
    def __init__(self):
        pass

    def __str__(self):
        s = f''
        return s
