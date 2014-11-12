

import types
import math
from math import cos, sin, sqrt
from vector3 import Vector3 as _vec3
from matrix44 import Matrix44 as _mat4

# Comparison threshold
_epsilon = 1E-9

class Quaternion(object):
    """Quaternion class.

    Quaternions are an extension to complex numbers and can be used
    to store rotations. They are composed of four floats which can be
    seen as an angle and an axis of rotation.
    """

    __slots__ = ('_q',)

    def __init__(self, *args):
        """Constructor.

        0 arguments: zeroes
        1 float argument:  w component, x,y,z = (0,0,0)
        4 arguments: components w,x,y,z

        If you want to construct from other types, use class methods.
        """

        # 0 arguments
        if not args:
            self._q = [ 0.0, 0.0, 0.0, 0.0 ]

        # 1 arguments
        elif len(args)==1:
            self._q = [ float(args[0]), 0.0, 0.0, 0.0 ]

        # 4 arguments
        elif len(args)==4:
            w,x,y,z = args
            self._q = map(float, args[:4])

        else:
            raise ValueError("Quaternion.__init__ takes 0, 1 or 4 parameters. ")


    @classmethod
    def blank(cls):
        q = cls.__new__(cls, object)
        q._q = [ 0.0, 0.0, 0.0, 0.0 ]
        return q

    @classmethod
    def identity(cls):
        q = cls.__new__(cls, object)
        q._q = [ 1.0, 0.0, 0.0, 0.0 ]
        return q

    @classmethod
    def from_angle_axis(cls, angle, axis):
        """Initialize self from an angle (in radians) and an axis and returns self."""

        q = cls.__new__(cls, object)

        if not any(axis):           # if the axis is (0.0, 0.0, 0.0)
            q._q = [ 1.0, 0.0, 0.0, 0.0 ]
            return q

        ax, ay, az = axis
        angle /= 2.0
        s = sin(angle) / sqrt(ax*ax + ay*ay + az*az)
        q._q = [ cos(angle), ax*s, ay*s, az*s ]
        return q

    @classmethod
    def from_floats(cls, w, x, y, z):
        """Creates a quaternion from individual float values.
        Warning: There is no checking (for efficiency) here: x, y, z _must_ be
        floats.
        """

        q = cls.__new__(cls, object)
        q._q = [w*1.0, x*1.0, y*1.0, z*1.0]
        return q

    @classmethod
    def from_two_vec3(cls, vstart, vend, dotprod=False, unit=False ):
        """ Returns NOT a UNIT quaternion! You have to normalise it yourself if you need to.
        Returns a quaternion describing a rotation from start to end:
            q.rotate( v_start )  = v_end

        'dotprod' is the dot product of two vectors. (optimisation if you calculate it beforehand)
        'unit' boolean, can be set to True if both vectors are of unit length.
        """
        q = cls.__new__(cls, object)

        if not dotprod:
            dotprod = vstart.dot(vend)

        if unit:
            w = 1.0 + dotprod
        else:
            w = sqrt(vstart.length_squared() * vend.length_squared()) + dotprod

        x, y, z = vstart.cross(vend)
        q._q = [ w, x, y, z ]
        return q

    @classmethod
    def from_euler_angles(cls, x, y, z):
        """ Create quaternion from 3 angles around axis given in radians. """

        q = cls.__new__(cls, object)

        x, y, z = x/2.0, y/2.0, z/2.0
        cx = cos(x)
        sx = sin(x)
        cy = cos(y)
        sy = sin(y)
        cz = cos(z)
        sz = sin(z)

        q._q = [ cy*cz*cx - sy*sz*sx,
                sy*sz*cx + cy*cz*sx,
                sy*cz*cx + cy*sz*sx,
                cy*sz*cx - sy*cz*sx ]
        return q

    @classmethod
    def from_matrix(cls, m):
        """Initialise self from either a mat3 or mat4 and returns self. """

        q = cls.__new__(cls, object)

        d1, d2, d3 = m[0,0], m[1,1], m[2,2]
        t = d1 + d2 + d3 + 1.0
        if t > _epsilon:            # most common case
            s = 0.5 / sqrt( t )
            q._q =[ 0.25/s,
                    (m[1,2] - m[2,1])*s,
                    (m[2,0] - m[0,2])*s,
                    (m[0,1] - m[1,0])*s ]
            return q

        if d1 >= d2 and d1 >= d3:
            s = sqrt(1.0 + d1 - d2 - d3)*2.0
            q._q =[ (m[2,1] + m[1,2])/s,
                    0.5/s,
                    (m[1,0] + m[0,1])/s,
                    (m[2,0] + m[0,2])/s ]
        elif d2 >= d1 and d2 >= d3:
            s = sqrt(1.0 + d2 - d1 - d3)*2.0
            q._q =[ (m[2,0] + m[0,2])/s,
                    (m[1,0] + m[0,1])/s,
                    0.5/s,
                    (m[2,1] + m[1,2])/s ]
        else:
            s = sqrt(1.0 + d3 - d1 - d2)*2.0
            q._q =[ (m[1,0] + m[0,1])/s,
                    (m[2,0] + m[0,2])/s,
                    (m[2,1] + m[1,2])/s,
                    0.5/s ]
        return q

    @classmethod
    def from_iter(cls, iterable):
        """Creates a quaternion from an iterable containing at least 4 values."""
        next = iter(iterable).next
        q = cls.__new__(cls, object)
        q._q = [ float(next()), float(next()), float(next()), float(next()) ]
        return q

    @classmethod
    def _from_float_sequence(cls, sequence):
        q = cls.__new__(cls, object)
        q._q = list(sequence[:4])
        return q

    def copy(self):
        """Returns an independent copy of this quaternion."""

        q = self.__new__(self.__class__, object)
        q._q = self._q[:]
        return q

    __copy__ = copy

    def _get_w(self):
        return self._q[0]
    def _set_w(self, w):
        try:
            self._q[0] = 1.0 * w
        except:
            raise TypeError, "Must be a number"
    w = property(_get_w, _set_w, None, "w component.")

    def _get_x(self):
        return self._q[1]
    def _set_x(self, x):
        try:
            self._q[1] = 1.0 * x
        except:
            raise TypeError, "Must be a number"
    x = property(_get_x, _set_x, None, "x component.")

    def _get_y(self):
        return self._q[2]
    def _set_y(self, y):
        try:
            self._q[2] = 1.0 * y
        except:
            raise TypeError, "Must be a number"
    y = property(_get_y, _set_y, None, "y component.")

    def _get_z(self):
        return self._q[3]
    def _set_z(self, z):
        try:
            self._q[3] = 1.0 * z
        except:
            raise TypeError, "Must be a number"
    z = property(_get_z, _set_z, None, "z component.")


    def _get_length(self):
        w, x, y, z = self._q
        return sqrt(w*w + x*x + y*y + z*z)

    def _set_length(self, length):
        q = self._q

        w, x, y, z = q
        l = length / sqrt(w*w + x*x + y*y +z*z)

        q[0] = w*l
        q[1] = x*l
        q[2] = y*l
        q[3] = z*l
        return self

    length = property(_get_length, _set_length, None, "Length of the quaternion")


    def set(self, w, x, y, z):
        """Sets the components of THIS quaternion.
        w -- w component
        x -- x component
        y -- y component
        z -- z component

        """
        q = self._q
        q[0] = w * 1.0
        q[1] = x * 1.0
        q[2] = y * 1.0
        q[3] = z * 1.0
        return self


    def __str__(self):

        w, x, y, z = self._q
        template = "({: f}, {: f}, {: f}, {: f})"
        return template.format( w, x, y, z )


    def __repr__(self):

        w, x, y, z = self._q
        return "Quaternion(%s, %s, %s, %s)" % (w, x, y, z)


    def __len__(self):

        return 4

    def __iter__(self):
        """Iterates the components in w, x, y, z order."""
        return iter(self._q[:])

    def __getitem__(self, index):
        """Retrieves a component, given its index.

        index -- 0, 1, 2 or 3 for w, x, y or z

        """
        try:
            return self._q[index]
        except IndexError:
            raise IndexError, "There are 4 values in this object, index should be 0, 1, 2 or 3!"

    def __setitem__(self, index, value):
        """Sets a component, given its index.

        index -- 0, 1, 2 or 3 for w, x, y or z
        value -- New (float) value of component

        """

        try:
            self._q[index] = 1.0 * value
        except IndexError:
            raise IndexError, "There are 4 values in this object, index should be 0, 1, 2 or 3!"
        except TypeError:
            raise TypeError, "Must be a number"


    def __eq__(self, rhs):
        """Test for equality
        rhs --quaternion sequence of 4 values
        """
        w, x, y, z = self._q
        ww, xx, yy, zz = rhs
        return abs(x-xx) < _epsilon and abs(y-yy) < _epsilon and abs(z-zz) < _epsilon


    def __ne__(self, rhs):
        """Test of inequality
        rhs --quaternionor sequenece of 4 values
        """
        x, y, z = self._q
        xx, yy, zz = rhs
        return abs(x-xx) > _epsilon or abs(y-yy) > _epsilon or abs(z-zz) > _epsilon


    def __hash__(self):

        return hash(self._q)

    def __add__(self, rhs):
        """Returns the result of adding a quaternion (or collection of 4 numbers)
        from this quaternion.

        rhs --quaternion or sequence of 4 values

        """

        w, x, y, z = self._q
        ow, ox, oy, oz = rhs
        return self.from_floats(w+ow, x+ox, y+oy, z+oz)


    def __iadd__(self, rhs):
        """Adds another quaternion (or a collection of 4 numbers) to this vector.

        rhs --quaternion or sequence of 4 values

        """
        ow, ox, oy, oz = rhs
        q = self._q
        q[0] += ow
        q[1] += ox
        q[2] += oy
        q[3] += oz
        return self


    def __sub__(self, rhs):
        """Returns the result of subtracting a quaternion (or collection of
        4 numbers) from this quaternion.

        rhs -- 4 values

        """

        w, x, y, z = self._q
        ow, ox, oy, oz = rhs
        return self.from_floats(w-ow, x-ox, y-oy, z-oz)


    def __isub__(self, rhs):
        """Subtracts another quaternion (or a collection of 4 numbers) from this
        quaternion.

        rhs --quaternion or sequence of 4 values

        """

        ow, ox, oy, oz = rhs
        q = self._q
        q[0] -= ow
        q[1] -= ox
        q[2] -= oy
        q[3] -= oz
        return self


    def __mul__(self, rhs):
        """Return the result of multiplying this quaternion by another quaternion
        or a sequence of 4 values.

        rhs -- Quaternion

        """
        w, x, y, z = self._q
        ow, ox, oy, oz = rhs
        return self.from_floats( w*ow - x*ox - y*oy - z*oz,
                                   w*ox + x*ow + y*oz - z*oy,
                                   w*oy + y*ow - x*oz + z*ox,
                                   w*oz + z*ow + x*oy - y*ox )


    def __imul__(self, rhs):
        """Multiply this quaternion by another quaternion or a sequence of 4 values.

        rhs -- Quaternion.

        """
        w, x, y, z = self._q
        ow, ox, oy, oz = rhs
        self._q = [ w*ow - x*ox - y*oy - z*oz,
                   w*ox + x*ow + y*oz - z*oy,
                   w*oy + y*ow - x*oz + z*ox,
                   w*oz + z*ow + x*oy - y*ox ]

        return self


    def __div__(self, rhs):
        """Return the result of dividing this quaternion by a scalar (single number)."""
        w, x, y, z = self._q
        return self.from_floats(w / rhs, x / rhs, y / rhs, z / rhs)


    def __idiv__(self, rhs):
        """Divide this quaternion by a scalar (single number)."""
        # maybe do: return self * other.get_inversed()
        q = self._q
        q[0] /= scalar
        q[0] /= scalar
        q[1] /= scalar
        q[2] /= scalar
        return self



    def __pow__(self, other):
        """Return self**quat."""

        q = other.get_scalar_multiplied( self.log() )
        return q.exp()

    def __neg__(self):
        """Negation. """

        w, x, y, z = self._q
        return self.from_floats(-w, -x, -y, -z)


    def __abs__(self):
        """Return magnitude.

        >>> q=Quaternion(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print round(abs(q),5)
        1.0
        """
        return self.length

     # a != 0
    def __nonzero__(self):
        return any( self._q )

    def scalar_mul(self, scalar):
        "Affects THIS quaternion"
        q = self._q
        q[0] *= scalar
        q[1] *= scalar
        q[2] *= scalar
        q[3] *= scalar

    def quat_mul(self, quat):
        "Affects THIS quaternion"
        w1, x1, y1, z1 = self._q
        w2, x2, y2, z2 = quat
        self._q = [ w1*w2 - x1*x2 - y1*y2 - z1*z2,
                   w1*x2 + x1*w2 + y1*z2 - z1*y2,
                   w1*y2 + y1*w2 - x1*z2 + z1*x2,
                   w1*z2 + z1*w2 + x1*y2 - y1*x2 ]
        return self

    def get_scalar_multiplied(self, scalar):
        "An independent quaternion object based on self."
        w, x, y, z = self._q
        return self.from_floats(w*scalar, x*scalar, y*scalar, z*scalar)

    def get_quat_multiplied(self, quat):
        "An independent quaternion object based on self."
        w1, x1, y1, z1 = self._q
        w2, x2, y2, z2 = quat
        return self.from_floats( w1*w2-x1*x2-y1*y2-z1*z2,
                                   w1*x2+x1*w2+y1*z2-z1*y2,
                                   w1*y2+y1*w2-x1*z2+z1*x2,
                                   w1*z2+z1*w2+x1*y2-y1*x2 )

    def scalar_div(self, scalar):
        "Affects THIS quat"
        q = self._q
        q[0] /= scalar
        q[0] /= scalar
        q[1] /= scalar
        q[2] /= scalar

    def sequence_div(self, seq):
        "Affects THIS quat"
        n0, n1, n2, n3 = seq
        q = self._q
        q[0] /= n0
        q[1] /= n1
        q[2] /= n2
        q[3] /= n3

    def get_scalar_divided(self, scalar):
        "An independent vector based on self."
        w, x, y, z = self._q
        return self.from_floats(w / scalar, x / scalar, y / scalar, z / scalar)

    def get_sequence_divided(self, seq):
        "An independent vector based on self."
        w, x, y, z = self._q
        n0, n1, n2, n3 = seq
        return self.from_floats(w/n0, x/n1, y/n2, z/n3)



    def get_conjugated(self):
        """Return conjugate.

        >>> q=Quaternion(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q.conjugate()
        (0.9689, -0.2160, -0.1080, -0.0540)
        """
        w, x, y, z = self._q
        return self.from_floats( w, -x, -y, -z )

    def get_normalized(self):
        """Return normalized quaternion.

        >>> q=Quaternion(0.9, 0.5, 0.2, 0.3)
        >>> q=q.normalize()
        >>> print q
        (0.8250, 0.4583, 0.1833, 0.2750)
        >>> print abs(q)
        1.0
        """
        nlen = 1.0/self.length
        w, x, y, z = self._q
        return self.from_floats( w*nlen, x*nlen, y*nlen, z*nlen )
    get_normalised = get_normalized

    def get_inversed(self):
        """Return inverse. Note that for unit quaternions
        inverse is the same as conjugate, and that conjugate is faster.

        >>> q=Quaternion(0.9, 0.5, 0.2, 0.3)
        >>> print q.inverse()
        (0.7563, -0.4202, -0.1681, -0.2521)
        """
        w, x, y, z = self._q
        len_2 = w*w + x*x + y*y + z*z
        return self.get_conjugated()/len_2

    def get_scaled(self, t):
        """
        Equivalent to nlerp'ing from identity quaternion to self.
        when t = 0 returns identity quaternion
        when t = 1.0 returns self
        """

        w, x, y, z = self._q
        nw = 1.0 + (w - 1.0)*t
        return self.from_floats( nw, x*t, y*t, z*t )

    def conjugate(self):
        """Return conjugate. Operates in place!

        >>> q=Quaternion(0.9689, 0.2160, 0.1080, 0.0540)
        >>> print q.conjugate()
        (0.9689, -0.2160, -0.1080, -0.0540)
        """
        w, x, y, z = self._q
        self._q = [ w, -x, -y, -z ]
        return self

    def normalize(self):
        """Return normalized quaternion. Operates in place!

        >>> q=Quaternion(0.9, 0.5, 0.2, 0.3)
        >>> q=q.normalize()
        >>> print q
        (0.8250, 0.4583, 0.1833, 0.2750)
        >>> print abs(q)
        1.0
        """
        nlen = 1.0 / self.length
        w, x, y, z = self._q
        self._q = [ w*nlen, x*nlen, y*nlen, z*nlen ]
        return self
    normalise = normalize

    def inverse(self):
        """Return inverse. Operates in place!

        >>> q=Quaternion(0.9, 0.5, 0.2, 0.3)
        >>> print q.inverse()
        (0.7563, -0.4202, -0.1681, -0.2521)
        """

        len_2 = w*w + x*x + y*y + z*z
        w, x, y, z = self._q
        self._q = [ w/len_2, -x/len_2, -y/len_2, -z/len_2 ]
        return self

    def scale(self, t):
        """
        Equivalent to nlerp'ing from identity quaternion to self. Operates in place!
        when t = 0 returns identity quaternion
        when t = 1.0 returns self
        """

        w, x, y, z = self._q
        nw = 1.0 + (w - 1.0)*t
        self._q = [ nw, x*t, y*t, z*t ]
        return self

    def to_angle_axis(self):
        """Return angle (in radians) and rotation axis.

        >>> q=Quaternion(0.9, 0.5, 0.2, 0.3)
        >>> angle, axis = q.to_angle_axis()
        >>> print round(angle,4)
        1.2011
        >>> print axis
        (0.8111, 0.3244, 0.4867)
        """

        w, x, y, z = self.get_normalised()._q

        # Clamp nself.w to protect against numerical inaccuracies.
        w = max(min(w, 1.0), -1.0)
        w = math.acos( w )
        s = sin( w )

        if s < _epsilon:
            return ( 0.0, _vec3(0.0, 0.0, 0.0) )
        return ( 2.0*w, _vec3(x/s, y/s, z/s) )

    def to_matrix33(self):
        """Return rotation matrix as mat3."""
        w, x, y, z = self._q
        xx = 2.0*x*x
        yy = 2.0*y*y
        zz = 2.0*z*z
        xy = 2.0*x*y
        zw = 2.0*z*w
        xz = 2.0*x*z
        yw = 2.0*y*w
        yz = 2.0*y*z
        xw = 2.0*x*w
        return _mat3(1.0-yy-zz, xy+zw, xz-yw,
                     xy-zw, 1.0-xx-zz, yz+xw,
                     xz+yw, yz-xw, 1.0-xx-yy)

    def to_matrix44(self):
        """Return rotation matrix as mat44."""
        w, x, y, z = self._q
        xx = 2.0*x*x
        yy = 2.0*y*y
        zz = 2.0*z*z
        xy = 2.0*x*y
        zw = 2.0*z*w
        xz = 2.0*x*z
        yw = 2.0*y*w
        yz = 2.0*y*z
        xw = 2.0*x*w
        return _mat4(1.0-yy-zz, xy+zw, xz-yw, 0.0,
                     xy-zw, 1.0-xx-zz, yz+xw, 0.0,
                     xz+yw, yz-xw, 1.0-xx-yy, 0.0,
                     0.0, 0.0, 0.0, 1.0)

    def to_log(self):
        """Return the natural logarithm of self."""

        q = cls.__new__(cls, object)
        w, x, y, z = self._q
        b = sqrt( x*x + y*y + z*z )
        if b <= _epsilon:
            if self.w <= _epsilon:
                raise ValueError( "math domain error" )
            q._q = [ math.log(w), 0.0, 0.0, 0.0 ]
            return q

        t = math.atan2( b, w )

        ct = cos(t)
        if abs(ct) <= _epsilon:
            raise ValueError( "math domain error" )

        r = w / ct
        if r <= _epsilon:
            raise ValueError( "math domain error" )

        f = t / b
        q._q = [ math.log(r), f*x, f*y, f*z ]
        return q


    def to_exp(self):
        """Return the exponential of self."""

        q = cls.__new__(cls, object)
        w, x, y, z = self._q

        b = sqrt(x*x + y*y + z*z)
        if b <= _epsilon:
            q._q = [ math.exp(w), 0.0, 0.0, 0.0 ]
        else:
            f = sin(b)/b
            q._q = [ math.exp(w)*cos(b), f*x, f*y, f*z ]
        return q


    def dot(self, other):
        """Return the dot product of self and b. """

        w, x, y, z = self._q
        ow, ox, oy, oz = other
        return w*ow + x*ox + y*oy + z*oz


    def rotate_vec3(self, seq):
        """Transforms a sequence of 3 items and returns the result as a vector.

        The quaternion must be a unit quaternion.
        This operation is equivalent to turning v into a Quaternion, computing
        self*v*self.conjugate() and turning the result back into a vec3.
        """
        vx, vy, vz = seq
        w, x, y, z = self._q
        ww = w*w
        xx = x*x
        yy = y*y
        zz = z*z
        wx = w*x
        wy = w*y
        wz = w*z
        xy = x*y
        xz = x*z
        yz = y*z

        return _vec3(ww*vx + xx*vx - yy*vx - zz*vx + 2*((xy-wz)*vy + (xz+wy)*vz),
                     ww*vy - xx*vy + yy*vy - zz*vy + 2*((xy+wz)*vx + (yz-wx)*vz),
                     ww*vz - xx*vz - yy*vz + zz*vz + 2*((xz-wy)*vx + (yz+wx)*vy))


    def rotate(self, seq):
        """Transforms a sequence of 3 values.

        The quaternion must be a unit quaternion.
        This operation is equivalent to turning v into a Quaternion, computing
        self*v*self.conjugate(), and returning the result as a tuple.
        """
        vx, vy, vz = seq
        w, x, y, z = self._q
        ww = w*w
        xx = x*x
        yy = y*y
        zz = z*z
        wx = w*x
        wy = w*y
        wz = w*z
        xy = x*y
        xz = x*z
        yz = y*z

        return tuple(ww*vx + xx*vx - yy*vx - zz*vx + 2*((xy-wz)*vy + (xz+wy)*vz),
                     ww*vy - xx*vy + yy*vy - zz*vy + 2*((xy+wz)*vx + (yz-wx)*vz),
                     ww*vz - xx*vz - yy*vz + zz*vz + 2*((xz-wy)*vx + (yz+wx)*vy))


def slerp(t, q0, q1, shortest=True):
    """Spherical linear interpolation between two quaternions.

    The return value is an interpolation between q0 and q1. For t=0.0
    the return value equals q0, for t=1.0 it equals q1.
    q0 and q1 must be unit quaternions.
    If shortest is True the interpolation is always done along the
    shortest path.
    """
    global _epsilon

    ca = q0.dot(q1)
    if shortest and ca<0:
        ca = -ca
        neg_q1 = True
    else:
        neg_q1 = False
    o = math.acos(ca)
    so = sin(o)

    if (abs(so)<=_epsilon):
        return Quaternion(q0)

    a = sin(o*(1.0-t)) / so
    b = sin(o*t) / so
    if neg_q1:
        return q0*a - q1*b
    else:
        return q0*a + q1*b


def squad(t, a, b, c, d):
    """Spherical cubic interpolation."""
    return slerp(2*t*(1.0-t), slerp(t,a,d), slerp(t,b,c))


def nlerp(t, q0, q1, shortest=True):
    """Normalised linear interpolation between two quaternions.

    The return value is an interpolation between q0 and q1. For t=0.0
    the return value equals q0, for t=1.0 it equals q1.
    q0 and q1 must be unit quaternions.
    If shortest is True the interpolation is always done along the
    shortest path.
    dp>0: short path
    dp<0: long path
    """
    dp = q0.dot(q1)
    if shortest and dp<0:
        q0 = -q0

    return Quaternion(q0.w + (q1.w-q0.w)*t,
                        q0.x + (q1.x-q0.x)*t,
                        q0.y + (q1.y-q0.y)*t,
                        q0.z + (q1.z-q0.z)*t).get_normalized()


def get_quat_to(v1, v2, dotprod=False, unit=False):
    """
    Returns a quaternion describing a rotation such that:
        q.rotate_vector(v1)  = v2

    So a rotation from v1 to v2.
    'dotprod' is the dot product of two vectors. (optimisation if you calculate it beforehand)
    'unit' boolean, are both vectors of unit length?
    """

    c = v1.cross(v2)
    if not any(c):
        return Quaternion(1.0, 0.0, 0.0, 0.0)

    x, y, z = c
    if not dotprod:
        dotprod = v1.dot(v2)

    if unit:
        w = 1 + dotprod
    else:
        w = sqrt(v1.length_squared() * v2.length_squared()) + dotprod

    return Quaternion(w, x, y, z).get_normalised()



def scale( t, quat ):
    """
    Equivalent to nlerp'ing from identity quaternion to 'quat'.
    when t = 0 returns identity quaternion
    when t = 1.0 returns 'quat'
    """

    w = 1.0 + (quat.w - 1.0)*t
    return Quaternion( w, quat.x*t, quat.y*t, quat.z*t )


