class HMM:
    '''A HMM class'''
    _pi = 0.0
    _A = 0.0
    _B = 0.0

    def __init__(self):
        self._pi = 5.0
        self._A = 0.0
        self._B = 0.0

    def get_pi(self):
        return self._pi

    def get_A(self):
        return self._A

    def get_B(self):
        return self._B
    
    def display(self):
        print 'pi =', self._pi
        print 'A =', self._A
        print 'B =', self._B
