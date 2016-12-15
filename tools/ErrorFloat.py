import math

def add_in_quad(a):
    return math.sqrt( sum([ x*x for x in a]) )

class ErrorFloat(object):
    def __init__(self, x=0.0, dx=None):
        if isinstance(x, ErrorFloat):
            ## copy constructor
            self.val = x.val
            self.errors = x.errors
        else:
            ## default constructor
            self.val = float(x)
            if dx is None:
                self.errors = [ math.sqrt(self.val) ]
            else:
                if isinstance(dx, list):
                    self.errors = dx
                else:
                    self.errors = [dx]

    def get_err(self):
        return add_in_quad(self.errors)

    err = property(get_err)

    def __repr__(self):
        strs = [ '%.4g' % x for x in [self.val] + self.errors ]
        return ' +/- '.join(strs)

    def latex(self):
        strs = [ '%.4g' % x for x in [self.val] + self.errors ]
        return  '$' + r' \pm '.join(strs) + '$'

    def __iadd__(self, other):
        assert isinstance(other, ErrorFloat)
        assert len(self.errors) == len(other.errors)
        self.errors = [ add_in_quad([e1, e2]) for e1, e2 in zip(self.errors, other.errors) ]
        self.val = self.val + other.val
        return self

    def __add__(self, other):
        if isinstance(other, ErrorFloat):
            x = ErrorFloat(self)
            x += other
            return x
        return self.val + other

    def __radd__(self, other):
        return self.__add__(other)

    def __isub__(self, other):
        assert isinstance(other, ErrorFloat)
        assert len(self.errors) == len(other.errors)
        self.errors = [ add_in_quad([e1, e2]) for e1, e2 in zip(self.errors, other.errors) ]
        self.val = self.val - other.val
        return self

    def __sub__(self, other):
        if isinstance(other, ErrorFloat):
            x = ErrorFloat(self)
            x -= other
            return x
        return self.val - other

    def __rsub__(self, other):
        return self.__sub__(other)

    def __imul__(self, other):
        assert isinstance(other, ErrorFloat)
        assert len(self.errors) == len(other.errors)
        self.errors = [ add_in_quad([e1*other.val, self.val*e2]) for e1, e2 in zip(self.errors, other.errors) ]
        self.val = self.val * other.val
        return self

    def __mul__(self, other):
        if isinstance(other, ErrorFloat):
            x = ErrorFloat(self)
            x *= other
            return x
        return self.val * other

    def __idiv__(self, other):
        assert isinstance(other, ErrorFloat)
        assert other.val, 'ErrorFloat.__idiv__ divide by zero.'
        assert len(self.errors) == len(other.errors)
        self.errors = [ add_in_quad([e1*other.val, self.val*e2]) / (other.val*other.val) for e1, e2 in zip(self.errors, other.errors) ]
        self.val = self.val / other.val
        return self

    def __div__(self, other):
        if isinstance(other, ErrorFloat):
            x = ErrorFloat(self)
            x /= other
            return x
        return self.val / other

    def add_err(self, e, index=0):
        len_errors = len(self.errors)
        if index >= len_errors:
            self.errors.extend( [0.0]*(index-len_errors+1) )
            print '  ErrorFloat.add_err: adding %s new error categories.' % (index-len_errors+1)
        self.errors[index] = add_in_quad([self.errors[index], e])
        return self

    def add_rel_err(self, e, index=0):
        return self.add_err(e*self.val, index)

    def add_sys(self, e):
        return self.add_err(e, 1)

    def add_rel_sys(self, e):
        return self.add_rel_err(e, 1)

    def make_stat_sys(self):
        len_errors = len(self.errors)
        assert len_errors >= 1
        if len_errors < 2:
            self.errors.append(0.0)
        self.errors[1] = add_in_quad([self.errors[0], self.errors[1]]) # sys = stat + sys
        self.errors[0] = 0.0 # stat = 0
        return self

