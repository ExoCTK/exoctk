# ! /usr/bin/env python
# quaternion module
"""Version 4 September 9, 2010 WMK
Flipped sign of the angle in the QX, QY, QZ, QJX, QJY, QJZ, set_values,
set_as_QX, ... functions to be consistent with the corrected multiplication.
Also updated
the doc strings.

Version 3 September 8, 2010 RLH
Backed out change to cnvrt in version 2.

Version 2 September 3, 2010 RLH
Fixed sign error in quaternion multiplication functions.
Added quaternion __str__ method.
Modified cnvrt method to return a CelestialVector


4/9/2010 WMK
Redefined the __init__ inputs
Changed from a Vector internal representation to 3 scalers
Fixed an error in cvt_att_Q_to_angles, was assuming an att2inertial Quaternion!
Streamlined some of the functions
 Version 1.0 August 3, 2010
 Got rid of degrees trig functions.

Combined this and rotationsx.py module to avoid circular imports and made it
PEP Compliant
Joe Filippazzo - 2018/06/26
"""
from math import radians, asin, cos, sin, sqrt, pi, degrees, atan2

from . import math_extensionsx as math2
from .astro_funcx import unit_limit

D2R = pi/180.
R2D = 180./pi
PI2 = 2.*pi


class GalacticPole:
    """Represents coordinates of galactic pole."""

    def __init__(self, latitude, longitude, ascending_node):
        """Initializes the coordinates of the galactic pole.

        Parameters
        ----------
        latitude: float
            Latitude of pole, in degrees.
        longitude: float
            Longitude of pole, in degrees.
        ascending_node: float
            Ascending node of pole, in degrees.
        """
        # Arguments specified in degrees, but values represented in radians.
        self.latitude = radians(latitude)
        self.longitude = radians(longitude)
        self.anode = radians(ascending_node)

    def __str__(self):
        """Returns string representation of the galactic pole."""
        text = (degrees(self.latitude), degrees(self.longitude),
                degrees(self.anode))
        # Convert attributes back into degrees for readability.
        return """GalacticPole: latitude: %.3fD, longitude: %.3fD,\
                  anode: %.3fD""" % text


# supports transformation to galactic coordinates
NGP = GalacticPole(192.859508, 27.128336, 32.932)


def QX(angle):
    """Creates rotation quaternion about X axis, rotates a vector about
    this axis.

    Parameters
    ----------
    angle: float
        The angle to rotate by.

    Result
    ------
    Quarternion
        The rotated quaternion.
    """
    return Quaternion(Vector(sin(angle/2.), 0., 0.), cos(angle/2.))


def QY(angle):
    """Creates rotation quaternion about Y axis, rotates a vector about
    this axis.

    Parameters
    ----------
    angle: float
        The angle to rotate by.

    Result
    ------
    Quarternion
        The rotated quaternion.
    """
    return Quaternion(Vector(0., sin(angle/2.), 0.), cos(angle/2.))


def QZ(angle):
    """Creates rotation quaternion about Z axis, rotates a vector about
    this axis

    Parameters
    ----------
    angle: float
        The angle to rotate by.

    Result
    ------
    Quarternion
        The rotated quaternion.
    """
    return Quaternion(Vector(0., 0., sin(angle/2.)), cos(angle/2.))


def Qmake_a_point(V):
    """Creates a pure Q, i.e. defines a pointing not a rotation

    Parameters
    ----------
    V: Vector
        The vector.

    Returns
    -------
    Quaternion
        The point as a quaternion.
    """
    return Quaternion(V, 0.)


def cvt_pt_Q_to_V(Q):
    """Converts a pure (pointing) Q to a unit position Vector

    Parameters
    ----------
    Q: Quaternion
        The quaternion to convert to a Vector.

    Returns
    -------
    Vector
        The point as a vector.
    """
    return Vector(Q.q1, Q.q2, Q.q3)


# The following functions are dependent upon the spacecraft definitions
# and perhaps should be moved to that module
def Qmake_body2inertial(coord1, coord2, V3pa):
    """Creates a rotation Q, going from the body frame to inertial.

    Parameters
    ----------
    coord1: float
        The first coordinate.
    coord2: float
        The second coordinate.
    V3pa: float
        The V3 position.

    Returns
    -------
    Quaternion
        The rotation quaternion
    """
    return QZ(coord1)*QY(-coord2)*QX(-V3pa)


def Qmake_v2v3_2body(v2, v3):
    """Creates a rotation Q, going from v2 and v3 in the body frame to
    inertial.

    Parameters
    ----------
    v2: float
        The V2 position.
    V3: float
        The V3 position.

    Returns
    -------
    Quaternion
        The rotation quaternion.
    """
    return QY(v3)*QZ(-v2)


def Qmake_v2v3_2inertial(coord1, coord2, V3pa, v2, v3):
    """Creates a rotation Q, going from v2 and v3 in the body frame to
    inertial

    Parameters
    ----------
    coord1: float
        The first coordinate.
    coord2: float
        The second coordinate.
    V3pa: float
        The V3 position.
    v2: float
        The V2 position.
    v3: float
        The V3 position.

    Returns
    -------
    Quaternion
        The rotation quaternion.
    """
    return QZ(coord1)*QY(-coord2)*QX(-V3pa)*QY(v3)*QZ(-v2)


def Qmake_aperture2inertial(coord1, coord2, APA, xoff, yoff, s, YapPA,
                            V3ref, V2ref):
    """Creates a rotation Q, going from the target in aperture frame to body.

    Parameters
    ----------
    coord1: float
        The first coordinate.
    coord2: float
        The second coordinate.
    APA: float
        The apature position.
    xoff: float
        The x offset.
    yoff: float
        The y offset.
    s: float
        The multiplicative factor.
    V2ref: float
        The V2 position.
    V3ref: float
        The V3 position.

    Returns
    -------
    Quaternion
        The rotation quaternion.
    """
    term1 = QZ(coord1)*QY(-coord2)*QX(-APA)*QY(-yoff)
    term2 = QZ(s*xoff)*QX(YapPA)*QY(V3ref)*QZ(-V2ref)
    return term1*term2


def cvt_body2inertial_Q_to_c1c2pa_tuple(Q):
    """Creates a angle tuple from Q, assuming body frame to inertial Q and
    321 rotation sequence.

    Parameters
    ----------
    Q: Quaternion
        The quaternion.

    Returns
    -------
    coord1 : float
        The first coordinate.
    coord2 : float
        The second coordinate.
    pa : float
        The poosition angle.
    """
    # Conversion from Euler symmetric parameters to matrix elements and
    # matrix elements to rotation angles is given in Isaac's papers
    r11 = Q.q1*Q.q1 - Q.q2*Q.q2 - Q.q3*Q.q3 + Q.q4*Q.q4
    r21 = 2.*(Q.q1*Q.q2 + Q.q3*Q.q4)
    r31 = 2.*(Q.q1*Q.q3 - Q.q2*Q.q4)
    r32 = 2.*(Q.q2*Q.q3 + Q.q1*Q.q4)
    r33 = -Q.q1*Q.q1 - Q.q2*Q.q2 + Q.q3*Q.q3 + Q.q4*Q.q4
    coord1 = atan2(r21, r11)
    if coord1 < 0.:
        coord1 += PI2
    coord2 = math2.asin2(r31)  # use "safe" version of sine
    pa = atan2(-r32, r33)
    if pa < 0.:
        pa += PI2
    return coord1, coord2, pa


def cvt_v2v3_using_body2inertial_Q_to_c1c2pa_tuple(Q, v2, v3):
    """Given Q and v2, v3 gives pos on sky and V3 PA.

    Parameters
    ----------
    Q: Quaternion
        The quaternion.
    v2: float
        The V2 position.
    v3: float
        The V3 position.

    Returns
    -------
    tuple
        The coordinates and position angle
    """
    Vp_body = Vector(0., 0., 0.)
    Vp_body.set_xyz_from_angs(v2, v3)
    Vp_eci_pt = Q.cnvrt(Vp_body)
    coord1 = atan2(Vp_eci_pt.y, Vp_eci_pt.x)
    if coord1 < 0.:
        coord1 += PI2
    coord2 = asin(unit_limit(Vp_eci_pt.z))

    V3_body = Vector(0., 0., 1.)
    V3_eci_pt = Q.cnvrt(V3_body)
    NP_eci = Vector(0., 0., 1.)
    V_left = cross(NP_eci, Vp_eci_pt)
    if V_left.length() > 0.:
        V_left = V_left/V_left.length()
    NP_in_plane = cross(Vp_eci_pt, V_left)
    x = dot(V3_eci_pt, NP_in_plane)
    y = dot(V3_eci_pt, V_left)
    pa = atan2(y, x)
    if pa < 0.:
        pa += PI2

    return coord1, coord2, pa


def cvt_c1c2_using_body2inertial_Q_to_v2v3pa_tuple(Q, coord1, coord2):
    """Given Q and a position, returns v2, v3, V3PA tuple

    Parameters
    ----------
    Q: Quaternion
        The quaternion.
    coord1: float
        The first coordinate.
    coord2: float
        The second coordinate

    Returns
    -------
    coord1 : float
        The first coordinate.
    coord2 : float
        The second coordinate.
    pa : float
        The poosition angle.
    """
    Vp_eci = Vector(1., 0., 0.)
    Vp_eci.set_xyz_from_angs(coord1, coord2)
    Vp_body_pt = Q.inv_cnvrt(Vp_eci)
    v2 = atan2(Vp_body_pt.y, Vp_body_pt.x)
    v3 = asin(unit_limit(Vp_body_pt.z))
    V3_body = Vector(0., 0., 1.)
    V3_eci_pt = Q.cnvrt(V3_body)
    NP_eci = Vector(0., 0., 1.)
    V_left = cross(NP_eci, Vp_eci)
    if V_left.length() > 0.:
        V_left = V_left / V_left.length()
    NP_in_plane = cross(Vp_eci, V_left)
    x = dot(V3_eci_pt, NP_in_plane)
    y = dot(V3_eci_pt, V_left)
    pa = atan2(y, x)
    if pa < 0.:
        pa += PI2
    return v2, v3, pa


class Quaternion:
    """This representation is used by Wertz and Markley"""
    def __init__(self, V, q4):
        """Quaternion constructor.

        Parameters
        ----------
        V: Vector
            The vector to construct the quaternion with.
        q4: Vector
            The fourth vector.
        """
        self.q1 = V.x
        self.q2 = V.y
        self.q3 = V.z
        self.q4 = q4

    def __str__(self):
        """Returns a string representation of the quaternion."""
        text = (self.q1, self.q2, self.q3, self.q4)
        return 'Quaternion: q1: %.3f, q2: %.3f, q3: %.3f, q4: %.3f' % text

    def length(self):
        """Returns length of the Q """
        a = self.q1*self.q1
        b = self.q2*self.q2
        c = self.q3*self.q3
        d = self.q4*self.q4
        return sqrt(a + b + c + d)

    def normalize(self):
        """Returns a copy of the Q normalized """
        scale = self.length()
        return Quaternion(Vector(self.q1/scale, self.q2/scale, self.q3/scale),
                          self.q4/scale)

    def conjugate(self):
        """Returns a copy of the conjugated Q """
        return Quaternion(Vector(-self.q1, -self.q2, -self.q3), self.q4)

    def __mul__(self, rs):
        """Defines Q*Q for quaternion multiplication.

        Parameters
        ----------
        rs: Quaternion
            The quaternion to multiply.

        Returns
        -------
        Q : Quaternion
            The multiplied Quaternion.
        """
        Q = Quaternion(Vector(0., 0., 0.), 0.)
        # Q.V = rs.V*self.q4 + self.V*rs.q4 + cross(self.V, rs.V)
        Q.q1 = rs.q1*self.q4 + self.q1*rs.q4 + (self.q2*rs.q3 - self.q3*rs.q2)
        Q.q2 = rs.q2*self.q4 + self.q2*rs.q4 + (self.q3*rs.q1 - self.q1*rs.q3)
        Q.q3 = rs.q3*self.q4 + self.q3*rs.q4 + (self.q1*rs.q2 - self.q2*rs.q1)
        Q.q4 = self.q4*rs.q4 - (self.q1*rs.q1 + self.q2*rs.q2 + self.q3*rs.q3)
        return Q

    def cnvrt(self, V):
        """Rotates a vector from the starting frame to the ending frame
        defined by the Q.

        Parameters
        ----------
        V: Vector
            The vector to rotate.

        Returns
        -------
        Vector
            The rotated Vector.
        """
        QV = Qmake_a_point(V)
        QV = self * QV * self.conjugate()
        return Vector(QV.q1, QV.q2, QV.q3)

    def inv_cnvrt(self, V):
        """Rotates a vector from the ending frame to the starting frame
        defined by the Q.

        Parameters
        ----------
        V: Vector
            The vector to invert.

        Returns
        -------
        Vector
            The inverted Vector.
        """
        QV = Qmake_a_point(V)
        QV = self.conjugate() * QV * self
        return Vector(QV.q1, QV.q2, QV.q3)

    def set_values(self, V, angle):
        """Sets quaterion values using a direction vector and a rotation of
        the coordinate frame about it.

        Parameters
        ----------
        V: Vector
            The direction Vector.
        angle: float
            The angle of rotation.
        """
        S = sin(-angle/2.)
        self.q1 = V.x * S
        self.q2 = V.y * S
        self.q3 = V.z * S
        self.q4 = cos(angle/2.)

    def set_as_QX(self, angle):
        """Sets quaterion in place like QX function.

        Parameters
        ----------
        angle: float
            The angle of rotation
        """
        self.q1 = sin(-angle/2.)
        self.q2 = 0.
        self.q3 = 0.
        self.q4 = cos(angle/2.)

    def set_as_QY(self, angle):
        """Sets quaterion in place like QY function

        Parameters
        ----------
        angle: float
            The angle of rotation
        """
        self.q1 = 0.
        self.q2 = sin(-angle/2.)
        self.q3 = 0.
        self.q4 = cos(angle/2.)

    def set_as_QZ(self, angle):
        """Sets quaterion in place like QZ function.

        Parameters
        ----------
        angle: float
            The angle of rotation.
        """
        self.q1 = 0.
        self.q2 = 0.
        self.q3 = sin(-angle/2.)
        self.q4 = cos(angle/2.)

    def set_as_mult(self, QQ1, QQ2):
        """Sets self as QQ1*QQ2 in place for quaternion multiplication.

        Parameters
        ----------
        QQ1: Quaternion
            The first quaternion.
        QQ2: Quaternion
            The second quaternion.
        """
        a = QQ1.q2*QQ2.q3 - QQ1.q3*QQ2.q2
        b = QQ1.q3*QQ2.q1 - QQ1.q1*QQ2.q3
        c = QQ1.q1*QQ2.q2 - QQ1.q2*QQ2.q1
        d = QQ1.q1*QQ2.q1 + QQ1.q2*QQ2.q2 + QQ1.q3*QQ2.q3
        self.q1 = QQ2.q1*QQ1.q4 + QQ1.q1*QQ2.q4 + a
        self.q2 = QQ2.q2*QQ1.q4 + QQ1.q2*QQ2.q4 + b
        self.q3 = QQ2.q3*QQ1.q4 + QQ1.q3*QQ2.q4 + c
        self.q4 = QQ1.q4*QQ2.q4 - d

    def set_as_point(self, V):
        """Set V as a point.

        Parameters
        ----------
        V: Vector
            The vector to set as a point.
        """
        self.q1 = V.x
        self.q2 = V.y
        self.q3 = V.z
        self.q4 = 0.

    def set_equal(self, Q):
        """Assigns values from other Q to this one.

        Parameters
        ----------
        Q: Quaternion
            The quaternion value to set.
        """
        self.q1 = Q.q1
        self.q2 = Q.q2
        self.q3 = Q.q3
        self.q4 = Q.q4

    def set_as_conjugate(self):
        """Assigns conjugate values in place. """
        self.q1 *= -1.
        self.q2 *= -1.
        self.q3 *= -1.


class NumericList(list):
    """List class that supports multiplication.  Only valid for numbers."""

    def __mul__(L1, L2):
        """Take the dot product of two numeric lists.
        Not using Vector for this because it is limited to three dimensions.
        Lists must have the same number of elements

        Parameters
        ----------
        L1: sequence
            The first list.
        L2: sequence
            The second list.

        Returns
        -------
        float
            The sum of the lists.
        """
        return(sum(map(lambda x, y: x*y, L1, L2)))


class Matrix(list):
    """Class to encapsulate matrix data and methods.

    A matrix is simply a list of lists that correspond to rows of the matrix.
    This is just intended to handle simple multiplication and vector rotations.
    For anything more advanced or computationally intensive, Python library
    routines should be used."""

    def __init__(self, rows):
        """Constructor for a matrix.

        This accepts a list of rows.
        It is assumed the rows are all of the same length.

        Parameters
        ----------
        rows: sequence
            The rows of the matrix.
        """
        for row in rows:
            self.append(NumericList(row))  # copy list

    def __str__(self):
        """Returns a string representation of the matrix."""
        return_str = 'Matrix:'

        for row_index in range(len(self)):
            row_str = 'Row %d: ' % (row_index + 1)
            row = self[row_index]

            for col_index in range(len(row)):
                row_str = row_str + '%6.3f  ' % (row[col_index])

            return_str = return_str + '\n' + row_str

        return(return_str)

    def element(self, row_index, col_index):
        """Returns an element of the matrix indexed by row and column.

        Indices begin with 0.

        Parameters
        ----------
        row_index: int
            The row index.
        col_index: int
            The column index.

        Returns
        -------
        float
            The matrix value.
        """
        return ((self[row_index])[col_index])

    def row(self, row_index):
        """Returns a specified row of the matrix.

        Parameters
        ----------
        row_index: int
            The row index.

        Returns
        -------
        list
            The row values.
        """

        return(self[row_index])

    def column(self, col_index):
        """Returns a specified column of the matrix as a numeric list.

        Parameters
        ----------
        col_index: int
            The column index.

        Returns
        -------
        list
            The column values.
        """
        return(NumericList([row[col_index] for row in self]))

    def num_rows(self):
        """Returns the number of rows in the matrix."""

        return(len(self))

    def num_cols(self):
        """Returns the number of columns in the matrix."""

        return (len(self[0]))  # assumes all rows of equal length

    def get_cols(self):
        """Returns list of all columns in a matrix."""
        rng = range(0, self.num_cols())
        return ([self.column(col_index) for col_index in rng])

    def __mul__(m1, m2):
        """Multiplies two Matrix objects and returns the resulting matrix.

        Number of rows in m1 must equal the number of columns in m2.

        Parameters
        ----------
        m1: Matrix
            The first matrix.
        m2: Matrix
            The second matrix.

        Returns
        -------
        Matrix
            The resultant matrix.
        """
        result_rows = []

        # Iterate over the rows in m1.  The first column of row i is formed by
        # multiplying the ith row of m1 by the first column of m2.  The second
        # column is formed by muliplying the ith row of m1 by the second
        # column of m2, etc.
        for row in m1:
            new_row = []

            for col in m2.get_cols():
                new_row.append(row * col)

            result_rows.append(new_row)

        return (Matrix(result_rows))


class Vector:
    "Class to encapsulate vector data and operations."

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Constructor for a three-dimensional vector.

        Note that two-dimensional vectors can be constructed by omitting one
        of the coordinates, which will default to 0.

        Parameters
        ----------
        x: float
            The x coordinate.
        y: float
            The y coordinate.
        z: float
            The z coordinate.
        """
        self.x = x     # Cartesian x coordinate
        self.y = y     # Cartesian y coordinate
        self.z = z     # Cartesian z coordinate

    def __str__(self):
        """Returns a string representation of the vector."""
        return('Vector: x: %.3f, y: %.3f, z: %.3f' % (self.x, self.y, self.z))

    def set_eq(self, x=None, y=None, z=None):
        """Assigns new value to vector.

        Arguments are now optional to permit this to be used with 2D vectors
        or to modify any subset of coordinates.

        Parameters
        ----------
        x: float
            The x coordinate.
        y: float
            The y coordinate.
        z: float
            The z coordinate.
        """
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if z is not None:
            self.z = z

    def length(self):
        """Returns magnitude of the vector """
        return(sqrt(self.x * self.x + self.y * self.y + self.z * self.z))

    def normalize(self):
        """Returns copy of the normalized vector """
        mag = self.length()
        return (Vector(self.x/mag, self.y/mag, self.z/mag))

    def __mul__(self, rs):
        """Implements Vector * scalar.  Can then use '*' syntax in multiplying
        a vector by a scalar rs

        Parameters
        ----------
        rs: float
            The scalar to multiply.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.x * rs
        y = self.y * rs
        z = self.z * rs
        return (Vector(x, y, z))

    def __rmul__(self, ls):
        """Implements float * Vector.

        Parameters
        ----------
        ls: float
            The scalar to multiply.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.x * ls
        y = self.y * ls
        z = self.z * ls
        return (Vector(x, y, z))

    def __add__(self, rs):
        """Implements Vector + Vector

        Parameters
        ----------
        rs: float
            The scalar to add.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.x + rs.x
        y = self.y + rs.y
        z = self.z + rs.z
        return (Vector(x, y, z))

    def __sub__(self, rs):
        """Implements Vector - Vector.

        Parameters
        ----------
        rs: float
            The scalar to subtract.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.x - rs.x
        y = self.y - rs.y
        z = self.z - rs.z
        return (Vector(x, y, z))

    def __truediv__(self, rs):
        """Implements Vector / float.

        Parameters
        ----------
        rs: float
            The scalar to divide.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.x / rs
        y = self.y / rs
        z = self.z / rs
        return (Vector(x, y, z))

    def __imul__(self, rs):
        """Implements Vector *= float

        Parameters
        ----------
        rs: float
            The scalar to multiply.

        Returns
        -------
        Vector
            The resultant vector.
        """
        self.x *= rs
        self.y *= rs
        self.z *= rs
        return (self)

    def __iadd__(self, rs):
        """Implements Vector += vector.

        Parameters
        ----------
        rs: float
            The scalar to add.

        Returns
        -------
        Vector
            The resultant vector.
        """
        self.x += rs.x
        self.y += rs.y
        self.z += rs.z
        return (self)

    def __isub__(self, rs):
        """Implements Vector -= vector.

        Parameters
        ----------
        rs: float
            The scalar to subtract.

        Returns
        -------
        Vector
            The resultant vector.
        """
        self.x -= rs.x
        self.y -= rs.y
        self.z -= rs.z
        return (self)

    def __idiv__(self, rs):
        """Implements Vector /= float.

        Parameters
        ----------
        rs: float
            The scalar to divide.

        Returns
        -------
        Vector
            The resultant vector.
        """
        self.x /= rs
        self.y /= rs
        self.z /= rs
        return (self)

    def create_matrix(self):
        """Converts a Vector into a single-column matrix."""

        column = [self.x, self.y, self.z]
        return(Matrix([[element] for element in column]))  # singleton list

    # Recommend deletion -- better to use a single interface that takes
    # two vectors.
    def dot(self, V2):
        """returns dot product between two vectors.

        Parameters
        ----------
        V2: Vector
            The vector to dot.

        Returns
        -------
        Vector
            The resultant vector.
        """
        return self.x * V2.x + self.y * V2.y + self.z * V2.z

    # Recommend deletion in favor of non-method version.
    def cross(self, V1, V2):
        """returns cross product of two vectors

        Parameters
        ----------
        V1: Vector
            The vector to cross.
        V2: Vector
            The vector to cross.

        Returns
        -------
        Vector
            The resultant vector.
        """
        x = self.y*V2.z - V1.z*V2.y
        y = self.z*V2.x - V1.x*V2.z
        z = self.x*V2.y - V1.y*V2.x
        return Vector(x, y, z)

    # Replace by separation - RLH
    def angle(self, V2):
        """Returns angle between the two vectors in degrees.

        Parameters
        ----------
        V2: Vector
            The vector to measure.

        Returns
        -------
        float
            The angle between the two vectors.
        """
        R1 = self.length()
        R2 = V2.length()
        adot = dot(self, V2)
        adot = adot / R1 / R2
        adot = min(1., adot)
        adot = max(-1., adot)
        return math2.acosd(adot)

    # RLH: What do these add?  We're creating methods just to access
    # individual attributes.
    def rx(self):
        """The magnitude of x"""
        return self.x

    def ry(self):
        """The magnitude of y"""
        return self.y

    def rz(self):
        """The magnitude of z"""
        return self.z

    # RLH: Suggest deletion in favor of __str__, which has the advantage
    # that it is called on print.
    def display(self):
        """Print the values"""
        return "[%f, %f, %f]" % (self.x, self.y, self.z)

    # RLH: Not necessary if CelestialVector is used.
    def set_xyz(self, ra, dec):
        """Creates a unit vector from spherical coordinates

        Parameters
        ----------
        ra: float
            The right ascension.
        dec: float
            The declination.
        """
        self.x = math2.cosd(dec) * math2.cosd(ra)
        self.y = math2.cosd(dec) * math2.sind(ra)
        self.z = math2.sind(dec)


class CelestialVector (Vector):
    "Class to encapsulate a unit vector on the celestial sphere."

    def __init__(self, ra=0.0, dec=0.0, frame='eq', degrees=True):
        """Constructor for a celestial vector.

        There are two spherical coordinates, a longitudinal coordinate (called
        right ascension), and a latitudinal coordinate (called declination).
        The RA is defined as the counterclockwise angle from a reference
        direction on the equatorial plane; it ranges from 0-360 degrees.
        The DEC is the angle between the vector and the equatorial plane;
        it ranges from -90 to 90 degrees. Angles are specified in degrees but
        represented internally as radians.

        The frame attribute indicates the coordinate frame of the vector,
        which may be 'eq' (equatorial, default), 'ec' (ecliptic), or 'gal'
        (galactic).  In equatorial coordinates, the equatorial plane is the
        celestial equator (extension of the Earth's equator) and the reference
        axis is the vernal equinox.  In ecliptic coordiantes, the equatorial
        plane is the ecliptic (the Earth's orbital plane) and the reference
        axis is usually defined relative to the Sun.  In galactic coordinates,
        the equatorial plane is the plane of the Galaxy.

        The degrees attribute should be True if the RA, DEC inputs are in
        degrees. Otherwise radians is assumed.

        The coordinates "ra" and "dec" may be used in all three systems.
        Other names for coordinates in different frames may be defined for
        clarity.

        A CelestialVector is also an ordinary unit vector, with Cartesian
        coordinates defined relative to the equatorial plane.

        Parameters
        ----------
        ra: float
            The right ascension.
        dec: float
            The declination.
        frame: str
            The frame to use.
        degrees: bool
            Use degrees.
        """
        if (degrees):
            ra = math2.D2R * ra
            dec = math2.D2R * dec

        self.ra = ra
        self.dec = dec
        self.frame = frame

        # Initialize standard vector with translated Cartesian coordinates
        x = cos(ra)*cos(dec)
        y = sin(ra)*cos(dec)
        z = sin(dec)
        Vector.__init__(self, x=x, y=y, z=z)

    def __str__(self, verbose=True):
        """Returns a string representation of the vector.  Displays angles
        in degrees.

        Parameters
        ----------
        verbose: bool
            Print some information.
        """
        a = (math2.R2D*self.ra, math2.R2D*self.dec, self.frame)
        celest_info = 'CelestialVector: RA: %.3fD, DEC: %.3fD, frame: %s' % a

        if (verbose):
            clss = super(CelestialVector, self).__str__()
            celest_info = celest_info + '\n' + clss
            return celest_info

    def set_eq(self, ra, dec, degrees=False):
        """Modifies a celestial vector with a new RA and DEC.

        degrees = True if units are degrees.  Default is radians.

        Parameters
        ----------
        ra: float
            The right ascension.
        dec: float
            The declination.
        degrees: bool
            Use degrees.
        """
        if (degrees):
            ra = math2.D2R * ra
            dec = math2.D2R * dec

            self.ra = ra
            self.dec = dec

        # Update Cartesian coordinates as well.
        x = cos(ra)*cos(dec)
        y = sin(ra)*cos(dec)
        z = sin(dec)
        super(CelestialVector, self).set_eq(x, y, z)

    def update_cartesian(self, x=None, y=None, z=None):
        """Modifies a celestial vector by specifying new Cartesian coordinates.

        Any subset of the Cartesian coordinates may be specifed.

        Parameters
        ----------
        x: float
            The extent in x.
        y: float
            The extent in y.
        z: float
            The extent in z.
        """

        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if z is not None:
            self.z = z

        self.ra = atan2(self.y, self.x)  # RA is arctan of y/x
        if (self.ra < 0):                 # Make sure RA is positive
            self.ra += 2*pi

            self.dec = math2.asin2(self.z)  # DEC is arcsin of z

    def rotate_about_axis(self, angle, axis):
        """This rotates a vector about an axis by the specified angle
        by using a rotation matrix.
        A new vector is returned.

        Axis must be 'x', 'y', or 'z'.
        The x-rotation rotates the y-axis toward the z-axis.
        The y-rotation rotates the z-axis toward the x-axis.
        The z-rotation rotates the x-axis toward the y-axis.

        Parameters
        ----------
        angle: float
            The angle of rotation.
        axis: str
            The axis to rotate about, ['x', 'y', 'z'].

        Returns
        -------
        result : Vector
            The rotated vector.
        """
        if (axis == 'x'):
            rot_matrix = Matrix([[1, 0, 0], [0, cos(angle), -sin(angle)],
                                 [0, sin(angle), cos(angle)]])

        elif (axis == 'y'):
            rot_matrix = Matrix([[cos(angle), 0, sin(angle)], [0, 1, 0],
                                 [-sin(angle), 0, cos(angle)]])

        elif (axis == 'z'):
            rot_matrix = Matrix([[cos(angle), -sin(angle), 0],
                                 [sin(angle), cos(angle), 0], [0, 0, 1]])

        else:
            print('Error')
            return

        new_matrix = rot_matrix * self.create_matrix()
        new_vector = new_matrix.column(0)
        result = CelestialVector()  # initialize with Cartesian coordiantes
        result.update_cartesian(x=new_vector[0], y=new_vector[1],
                                z=new_vector[2])
        return result

    def rotate_about_eigenaxis(self, angle, eigenaxis):
        """Rotates a vector about arbitrary eigenaxis.

        eigenaxis = Vector object (axis about which to rotate).
        angle = angle to rotate by in radians.
        Rotation is counterclockwise looking outward from origin along
        eigenaxis. Function uses rotation matrix from Rodrigues formula.

        Note: This function is more general than rotate_about_axis above and
        could be used in its place.  However, rotate_about_axis is faster and
        clearer when the rotation axis is one of the Cartesian axes.

        Parameters
        ----------
        angle: float
            The angle of rotation.
        eigenaxis: Vector
            The eigenaxis to rotate about.

        Returns
        -------
        result : Vector
            The rotated vector.
        """
        cos_ang = cos(angle)    # Used repeatedly below
        sin_ang = sin(angle)

        # Fill out the Rodrigues rotation matrix
        R11 = cos_ang + eigenaxis.x**2 * (1 - cos_ang)
        R12 = eigenaxis.x * eigenaxis.y * (1 - cos_ang) - eigenaxis.z * sin_ang
        R13 = eigenaxis.x * eigenaxis.z * (1 - cos_ang) + eigenaxis.y * sin_ang
        R21 = eigenaxis.x * eigenaxis.y * (1 - cos_ang) + eigenaxis.z * sin_ang
        R22 = cos_ang + eigenaxis.y**2 * (1 - cos_ang)
        R23 = eigenaxis.y * eigenaxis.z * (1 - cos_ang) - eigenaxis.x * sin_ang
        R31 = eigenaxis.x * eigenaxis.z * (1 - cos_ang) - eigenaxis.y * sin_ang
        R32 = eigenaxis.y * eigenaxis.z * (1 - cos_ang) + eigenaxis.x * sin_ang
        R33 = cos_ang + eigenaxis.z**2 * (1 - cos_ang)

        r1, r2, r3 = [R11, R12, R13], [R21, R22, R23], [R31, R32, R33]
        rot_matrix = Matrix([r1, r2, r3])
        new_matrix = rot_matrix * self.create_matrix()
        new_vector = new_matrix.column(0)
        result = CelestialVector()  # initialize with Cartesian coordinates
        result.update_cartesian(x=new_vector[0], y=new_vector[1],
                                z=new_vector[2])
        return(result)

    def rotate_using_quaternion(self, angle, eigenaxis):
        """Rotates a vector about arbitrary eigenaxis using quaternion.

        This is an alternative formulation for rotate_about_eigenaxis.
        Interface is the same as rotate_about_eigenaxis.

        Parameters
        ----------
        angle: float
            The angle of rotation.
        eigenaxis: Vector
            The eigenaxis to rotate about.

        Returns
        -------
        Vector
            The rotated vector.
        """
        q = Quaternion(eigenaxis, 0.0)

        # Need to negate here because set_values performs a negative rotation
        # quaternion now represents the rotation
        q.set_values(eigenaxis, -angle)
        return(make_celestial_vector(q.cnvrt(self)))

    def transform_frame(self, new_frame):
        """Transforms coordinates between celestial and ecliptic frames
        and returns result as a new CelestialVector.
        If new coordinate frame is the same as the old, a copy of the vector
        is returned.

        Parameters
        ----------
        new_frame: str
            Convert to new frame.

        Returns
        -------
        result : Vector
            The transformed vector.
        """
        result = None
        gal_ec = new_frame == 'gal' and self.frame == 'ec'
        ec_gal = new_frame == 'ec' and self.frame == 'gal'

        # Equatorial to ecliptic: rotate z-axis toward y-axis.
        if ((new_frame == 'ec') and (self.frame == 'eq')):
            result = self.rotate_about_axis(-math2.OBLIQUITY, 'x')

        # Ecliptic to equatorial: rotate y-axis toward z-axis.
        elif ((new_frame == 'eq') and (self.frame == 'ec')):
            result = self.rotate_about_axis(math2.OBLIQUITY, 'x')

        elif ((new_frame == 'gal') and (self.frame == 'eq')):
            # Use formula from Wayne Kinzel's book, adjusted for
            # J2000 coordinates.
            a = cos(self.dec)*cos(NGP.longitude)*cos(self.ra - NGP.latitude)
            b = math2.asin2(a + sin(self.dec) * sin(NGP.longitude))
            arg1 = sin(self.dec) - sin(b)*sin(NGP.longitude)
            arg2 = cos(self.dec)*sin(self.ra - NGP.latitude)*cos(NGP.longitude)

            lng = atan2(arg1, arg2) + NGP.anode

            result = CelestialVector(lng, b, degrees=False)

        elif ((new_frame == 'eq') and (self.frame == 'gal')):
            lng = self.ra   # use l, b notation here for clarity
            b = self.dec
            term1 = cos(b) * cos(NGP.longitude) * sin(lng - NGP.anode)
            dec = math2.asin2(term1 + sin(b) * sin(NGP.longitude))
            arg1 = cos(b) * cos(lng - NGP.anode)
            sinterm = sin(NGP.longitude) * sin(lng - NGP.anode)
            arg2 = sin(b) * cos(NGP.longitude) - cos(b) * sinterm
            ra = atan2(arg1, arg2) + NGP.latitude

            result = CelestialVector(ra, dec, degrees=False)

        elif gal_ec or ec_gal:
            print("""Error: Direct conversion between ecliptic and\
                     galactic coordinates not supported yet""")

        elif new_frame != self.frame:
            print("Error: unrecognized coordinate frame.")

        # If there was an error, return a copy of the initial vector.
        if result is None:
            result = CelestialVector(self.ra, self.dec, self.frame, False)

        else:
            result.frame = new_frame  # record new frame

        return (result)

    def rotate_by_posang(self, pa):
        """Returns the vector that results from rotating the self vector
        counterclockwise from the North projection onto the plane
        orthogonal to that vector by the specified position angle
        (in radians). See "V3-axis Position Angle", John Isaacs, May 2003 for
        further discussion.

        Parameters
        ----------
        pa: float
            The position angle.

        Returns
        -------
        result : Vector
            The rotated vector.
        """
        x_coord = -cos(self.ra)*sin(self.dec)*cos(pa) - sin(self.ra)*sin(pa)
        y_coord = -sin(self.ra)*sin(self.dec)*cos(pa) + cos(self.ra)*sin(pa)
        z_coord = cos(self.dec)*cos(pa)
        result = CelestialVector()
        result.update_cartesian(x_coord, y_coord, z_coord)
        return(result)

    def position_angle(self, v):
        """Returns the position angle of v at the self vector, in radians.

        v is an arbitrary vector that should be a CelestialVector object.
        The position angle is the angle between the North vector on the
        plane orthogonal to the self vector and the projection of v onto
        that plane, defined counterclockwise.
        See "V3-axis Position Angle", John Isaacs, May 2003 for
        further discussion.

        Parameters
        ----------
        v: Vector
            The vector to measure against.

        Returns
        -------
        pa : float
            The position angle between the two vectors.
        """
        y_coord = cos(v.dec) * sin(v.ra - self.ra)
        b = cos(v.dec) * sin(self.dec) * cos(v.ra - self.ra)
        x_coord = sin(v.dec) * cos(self.dec) - b
        pa = atan2(y_coord, x_coord)

        if (pa < 0):
            pa += (2 * pi)  # PA has range 0-360 degrees

        return(pa)


class Attitude(CelestialVector):
    "Defines an Observatory attitude by adding a position angle."""

    def __init__(self, ra=0.0, dec=0.0, pa=0.0, frame='eq', degrees=True):
        """Constructor for an Attitude.

        pa = position_angle in degrees(default) or radians if degrees=False
        is specified. Other arguments are the same as with CelestialVector

        Parameters
        ----------
        ra: float
            The right ascension.
        dec: float
            The declination.
        pa: float
            The position angle.
        frame: str
            The frame to use.
        degrees: bool
            Use degrees.
        """
        super(Attitude, self).__init__(ra=ra, dec=dec, frame=frame,
                                       degrees=degrees)

        if (degrees):   # convert into radians
            pa = math2.D2R * pa

        self.pa = pa

    def __str__(self, verbose=True):
        """Returns a string representation of the attitude.

        verbose (optional) = flag indicating whether detailed Vector
        information should be included.

        Parameters
        ----------
        verbose: bool
            Print information.
        
        Returns
        -------
        att_info : str
            A string representation of the attitute. 
        """
        att_info = 'Attitude: PA: %.3fD' % (math2.R2D * self.pa)
        att_info = att_info + '\n' + super(Attitude, self).__str__(verbose)
        return att_info


# Functions that operate on vectors but are not methods.
def dot(v1, v2):
    """returns dot product between two vectors, non class member.

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.

    Returns
    -------
    float
        The dot product of the vectors.
    """
    return(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)


def cross(v1, v2):
    """Returns cross product between two vectors, non class member.

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.

    Returns
    -------
    float
        The cross product of the vectors.
    """
    x = v1.y*v2.z - v1.z*v2.y
    y = v1.z*v2.x - v1.x*v2.z
    z = v1.x*v2.y - v1.y*v2.x
    return Vector(x, y, z)


def separation(v1, v2, norm=False):
    """Returns angle between two unit vectors in radians.

    The angle between two normalized vectors is the arc-cosine of the dot
    product. Unless the norm attribute is set to True, it is assumed the
    vectors are already normalized (for performance).

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.
    norm: bool
        Normalize the vectors.

    Returns
    -------
    separation : float
        The separation of the vectors.
    """
    if (norm):
        v1 = v1.normalize()
        v2 = v2.normalize()

    separation = math2.acos2(dot(v1, v2))

    # For very small angles, cos and acos behave poorly as the cosine of a
    # very small angle is interpreted as 1.0.  Therefore, recompute using
    # the cross product if the result is less than 1 degree.
    if (separation < math2.D2R):
        vcross = cross(v1, v2)
        separation = math2.asin2(vcross.length())

    return(separation)


def ra_delta(v1, v2):
    """Returns difference in right ascension between two CelestialVectors.

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.

    Returns
    -------
    delta_ra : float
        The difference in RA.
    """
    delta_ra = v1.ra - v2.ra

    # Check for zero crossings.  If the difference exceeds 180 degrees,
    # adjust by 360 in opposite direction.

    if (delta_ra < -pi):
        delta_ra = delta_ra + 2*pi
    elif (delta_ra > pi):
        delta_ra = delta_ra - 2*pi

    return(delta_ra)


def ra_separation(v1, v2):
    """Returns separation in right ascension between two CelestialVectors.
    This is accurate only if the difference in declination is small.

    |sep| = DELTA-RA cos DEC

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.

    Returns
    -------
    float
        The separation between RA values.
    """
    delta_ra = ra_delta(v1, v2)
    dec = math2.avg2(v1.dec, v2.dec)  # use average of the two declinations.
    return(delta_ra * cos(dec))


def dec_separation(v1, v2):
    """Returns difference in declination between two CelestialVectors.

    Parameters
    ----------
    v1: Vector
        The first vector.
    v2: Vector
        The second vector.

    Returns
    -------
    float
        The separation between Dec values
    """
    return(v1.dec - v2.dec)    # simply take the difference in declination


def make_celestial_vector(v):
    """Takes a Vector object and creates an equivalent CelestialVector.

    Input vector v must be a unit vector.

    Parameters
    ----------
    v: Vector
        The vector to convert.

    Returns
    -------
    result : Vector
        The updated vector.
    """
    result = CelestialVector()
    result.update_cartesian(v.x, v.y, v.z)
    return(result)


def projection(v, axis):
    """Returns projection of vector v on plane normal to axis.

    First take cross-product of v and the axis and normalize it.
    Then cross the axis with the result and return a CelestialVector.
    See http://www.euclideanspace.com/maths/geometry/elements/plane/
    lineOnPlane/index.htm.

    Parameters
    ----------
    v: Vector
        The vector to convert.
    axis: str
        The axis to project onto.

    Returns
    -------
    Vector
        The updated vector.
    """
    return(make_celestial_vector(cross(axis, (cross(v, axis)).normalize())))


def pos_V_to_ra_dec(V):
    """Returns tuple of spherical angles from unit direction Vector.

    Parameters
    ----------
    V: Vector
        The vector to analyze.

    Returns
    -------
    ra : float
        The ra of the vector.
    dec : float
        The dec of the vector.
    """
    ra = math2.atan2d(V.y, V.x)
    V.z = min(1., V.z)
    V.z = max(-1., V.z)
    dec = math2.asind(V.z)
    if ra < 0.:
        ra += 360.
    return(ra, dec)


# RLH: Recommend replacement by separation.
def angle(V1, V2):
    """returns angle between two vectors in degrees, non class member.

    Parameters
    ----------
    V1: Vector
        The first vector.
    V2: Vector
        The second vector.

    Returns
    -------
    float
        The angle between the vectors.
    """
    R1 = V1.length()
    R2 = V2.length()
    adot = dot(V1, V2)
    adot = adot / R1 / R2
    adot = min(1., adot)
    adot = max(-1., adot)
    return math2.acosd(adot)


def vel_ab(U, Vel):
    """Takes a unit vector and a velocity vector(km/s) and returns a unit
    vector modidifed by the velocity abberation.

    Parameters
    ----------
    U: Vector
        The unit vector.
    Vel: Vector
        The velocity vector to multiply.

    Returns
    -------
    Vector
        The modified vector.
    """
    c = 2.9979e5  # speed of light in km/s
    Beta = Vel * (1./c)
    rgamma = sqrt(1.-dot(Beta, Beta))  # This is 1/gamma
    ubeta = dot(U, Beta)
    b = (1./(1. + ubeta))
    return (U*rgamma + Beta * (1. + (1.-rgamma)*ubeta/dot(Beta, Beta))) * b
