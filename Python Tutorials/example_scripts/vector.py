
class Vector:
    
    def __init__(self, x=0, y=0, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __add__(self, v):
        x, y, z = self.x + v.x, self.y + v.y, self.z + v.z
        return self.__class__(x, y, z)        

    def __sub__(self, v):
        x, y, z = self.x - v.x, self.y - v.y, self.z - v.z
        return self.__class__(x, y, z)

    def __mul__(self, v):
        # does a dot product
        return self.x * v.x + self.y * v.y + self.z * v.z

    def __iadd__(self, v):
        self.x += v.x
        self.y += v.y
        self.z += v.z
    
    def __isub__(self, v):
        self.x -= v.x
        self.y -= v.y
        self.z -= v.z

    def __eq__(self, v):
        return self.magnitude == v.magnitude

    def __ge__(self, v):
        return self.magnitude >= v.magnitude
        
    def __gt__(self, v):
        return self.magnitude > v.magnitude
    
    def __le__(self, v):
        return self.magnitude <= v.magnitude        
        
    def __lt__(self, v):
        return self.magnitude < v.magnitude      
        
    def __ne__(self, v):
        return not self.magnitude == v.magnitude        
        
    def __repr__(self):
        return '%s at %i (%.3f, %.3f5, %.3f)' % (self.__class__.__name__, id(self), self.x, self.y, self.z)
        
    def __getitem__(self, i):
        if i > 2:
            raise IndexError('Vector is only 3D!')
        return (self.x, self.y, self.z)[i]

    @property
    def magnitude(self):
        return (self.x**2 + self.y**2 + self.z**2)**(0.5)
