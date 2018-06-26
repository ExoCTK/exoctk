#! /usr/bin/env python
#quaternion module

"""Version 4 September 9, 2010 WMK
Flipped sign of the angle in the QX, QY, QZ, QJX, QJY, QJZ, set_values, set_as_QX, ...
functions to be consistent with the corrected multiplication. Also updated
the doc strings."""

"""Version 3 September 8, 2010 RLH
Backed out change to cnvrt in version 2."""

"""Version 2 September 3, 2010 RLH
Fixed sign error in quaternion multiplication functions.
Added quaternion __str__ method.
Modified cnvrt method to return a CelestialVector."""


"""4/9/2010 WMK
Redefined the __init__ inputs
Changed from a Vector internal representation to 3 scalers
Fixed an error in cvt_att_Q_to_angles, was assuming an att2inertial Quaternion!
Streamlined some of the functions
 Version 1.0 August 3, 2010
 Got rid of degrees trig functions."""

"""
Combined this and rotationsx.py module to avoid circular imports and made it PEP Compliant
Joe Filippazzo - 2018/06/26
"""

import math
from . import math_extensionsx as math2

D2R = math.pi/180.
R2D = 180. / math.pi
PI2 = 2. * math.pi
unit_limit = lambda x: min(max(-1., x), 1.)

class GalacticPole (object):
    """Represents coordinates of galactic pole."""

    def __init__ (self, latitude, longitude, ascending_node):
        """Initializes the coordinates of the galactic pole.

        latitude = latitude of pole, in degrees
        longitude = longitude of pole, in degrees
        ascending_node = ascending node of pole, in degrees."""

        #Arguments specified in degrees, but values represented in radians.
        self.latitude = math.radians(latitude)
        self.longitude = math.radians(longitude)
        self.anode = math.radians(ascending_node)

    def __str__ (self):
        """Returns string representation of the galactic pole."""

        #Convert attributes back into degrees for readability.
        return('GalacticPole: latitude: %.3fD, longitude: %.3fD, anode: %.3fD'\
        %(degrees(self.latitude), degrees(self.longitude), degrees(self.anode)))

NGP = GalacticPole(192.859508, 27.128336, 32.932)
# supports transformation to galactic coordinates


def QX(angle):
   """Creates rotation quaternion about X axis, rotates a vector about this axis"""
   return Quaternion(Vector(math.sin(angle/2.), 0., 0.), math.cos(angle/2.))

def QY(angle):
   """Creates rotation quaternion about Y axis, rotates a vector about this axis"""
   return Quaternion(Vector(0., math.sin(angle/2.), 0.), math.cos(angle/2.))

def QZ(angle):
   """Creates rotation quaternion about Z axis, rotates a vector about this axis"""
   return Quaternion(Vector(0., 0., math.sin(angle/2.)), math.cos(angle/2.))

def Qmake_a_point(V):
   """Creates a pure Q, i.e. defines a pointing not a rotation"""
   return Quaternion(V, 0.)

def cvt_pt_Q_to_V(Q):
   """Converts a pure (pointing) Q to a unit position Vector"""
   return Vector(Q.q1, Q.q2, Q.q3)


#The following functions are dependent upon the spacecraft definitions and perhaps should be moved to that module

def Qmake_body2inertial(coord1, coord2, V3pa):
   """Creates a rotation Q, going from the body frame to inertial"""
   return QZ(coord1)*QY(-coord2)*QX(-V3pa)

def Qmake_v2v3_2body(v2, v3):
   """Creates a rotation Q, going from v2 and v3 in the body frame to inertial"""
   return QY(v3)*QZ(-v2)

def Qmake_v2v3_2inertial(coord1, coord2, V3pa, v2, v3):
   """Creates a rotation Q, going from v2 and v3 in the body frame to inertial"""
   return QZ(coord1)*QY(-coord2)*QX(-V3pa)*QY(v3)*QZ(-v2)

def Qmake_aperture2inertial(coord1, coord2, APA, xoff, yoff, s, YapPA, V3ref, V2ref):
   """Creates a rotation Q, going from the target in aperture frame to body"""
   return QZ(coord1)*QY(-coord2)*QX(-APA)*QY(-yoff)*QZ(s*xoff)*QX(YapPA)*QY(V3ref)*QZ(-V2ref)

def cvt_body2inertial_Q_to_c1c2pa_tuple(Q):
   """Creates a angle tuple from Q, assuming body frame to inertial Q and 321 rotation sequence"""
   #Conversion from Euler symmetric parameters to matrix elements and matrix elements to rotation angles
   # is given in Isaac's papers
   r11 = Q.q1*Q.q1 - Q.q2*Q.q2 - Q.q3*Q.q3 + Q.q4*Q.q4
   r21 = 2.*(Q.q1*Q.q2 + Q.q3*Q.q4)
   r31 = 2.*(Q.q1*Q.q3 - Q.q2*Q.q4)
   r32 = 2.*(Q.q2*Q.q3 + Q.q1*Q.q4)
   r33 = -Q.q1*Q.q1 - Q.q2*Q.q2 + Q.q3*Q.q3 + Q.q4*Q.q4
   coord1 = math.atan2(r21, r11)
   if coord1 < 0. : coord1 += PI2
   coord2 = math.asin2(r31)  # use "safe" version of sine
   pa  = math.atan2(-r32, r33)
   if pa < 0. : pa += PI2
   return (coord1, coord2, pa)

def cvt_v2v3_using_body2inertial_Q_to_c1c2pa_tuple(Q, v2, v3):
   """Given Q and v2, v3 gives pos on sky and V3 PA """
   Vp_body = Vector(0., 0., 0.)
   Vp_body.set_xyz_from_angs(v2, v3)
   Vp_eci_pt = Q.cnvrt(Vp_body)
   coord1  = math.atan2(Vp_eci_pt.y, Vp_eci_pt.x)
   if coord1 < 0. : coord1 += PI2
   coord2 = math.asin(unit_limit(Vp_eci_pt.z))

   V3_body = Vector(0., 0., 1.)
   V3_eci_pt = Q.cnvrt(V3_body)
   NP_eci = Vector(0., 0., 1.)
   V_left = math.cross(NP_eci, Vp_eci_pt)
   if V_left.length()>0.:
      V_left = V_left/V_left.length()
   NP_in_plane = math.cross(Vp_eci_pt, V_left)
   x = math.dot(V3_eci_pt, NP_in_plane)
   y = math.dot(V3_eci_pt, V_left)
   pa  = math.atan2(y, x)
   if pa < 0. : pa += PI2

   return (coord1, coord2, pa)

def cvt_c1c2_using_body2inertial_Q_to_v2v3pa_tuple(Q, coord1, coord2):
   """Given Q and a position, returns v2, v3, V3PA tuple """
   Vp_eci = Vector(1., 0., 0.)
   Vp_eci.set_xyz_from_angs(coord1, coord2)
   Vp_body_pt = Q.inv_cnvrt(Vp_eci)
   v2 = math.atan2(Vp_body_pt.y, Vp_body_pt.x)
   v3 = math.asin(unit_limit(Vp_body_pt.z))
   V3_body = Vector(0., 0., 1.)
   V3_eci_pt = self.cnvrt(V3_body)
   NP_eci = Vector(0., 0., 1.)
   V_left = math.cross(NP_eci, Vp_eci)
   if V_left.length()>0.:
    V_left = V_left / V_left.length()
   NP_in_plane = math.cross(Vp_eci, V_left)
   x = math.dot(V3_eci_pt, NP_in_plane)
   y = math.dot(V3_eci_pt, V_left)
   pa  = math.atan2(y, x)
   if pa < 0. : pa += PI2
   return (v2, v3, pa)


############################################################


class Quaternion:
   """This representation is used by Wertz and Markley """
   def __init__(self, V, q4):
      """Quaternion constructor """
      self.q1 = V.x
      self.q2 = V.y
      self.q3 = V.z
      self.q4 = q4

   def __str__(self):
       """Returns a string representation of the quaternion."""

       return('Quaternion: q1: %.3f, q2: %.3f, q3: %.3f, q4: %.3f'\
               % (self.q1, self.q2, self.q3, self.q4))

   def length(self):
      """Returns length of the Q """
      return math.math.sqrt(self.q1*self.q1 + self.q2*self.q2 + self.q3*self.q3 + self.q4*self.q4)

   def normalize(self):
      """Returns a copy of the Q normalized """
      scale = self.length()
      return Quaternion(Vector(self.q1/scale, self.q2/scale, self.q3/scale), self.q4/scale)

   def conjugate(self):
      """Returns a copy of the conjugated Q """
      return Quaternion(Vector(-self.q1, -self.q2, -self.q3), self.q4)

   def __mul__(self, rs):
      """Defines Q*Q for quaternion multiplication """
      Q = Quaternion(Vector(0., 0., 0.), 0.)
      #Q.V = rs.V * self.q4 + self.V * rs.q4 + cross(self.V, rs.V)
      Q.q1 = rs.q1 * self.q4 + self.q1 * rs.q4 + (self.q2 * rs.q3 - self.q3 * rs.q2)
      Q.q2 = rs.q2 * self.q4 + self.q2 * rs.q4 + (self.q3 * rs.q1 - self.q1 * rs.q3)
      Q.q3 = rs.q3 * self.q4 + self.q3 * rs.q4 + (self.q1 * rs.q2 - self.q2 * rs.q1)
      Q.q4 = self.q4 * rs.q4 - (self.q1 * rs.q1 + self.q2 * rs.q2 + self.q3 * rs.q3)
      return Q

   def cnvrt(self, V):
      """Rotates a vector from the starting frame to the ending frame defined by the Q """
      QV = Qmake_a_point(V)
      QV = self * QV * self.conjugate()
      return Vector(QV.q1, QV.q2, QV.q3)

   def inv_cnvrt(self, V):
      """Rotates a vector from the ending frame to the starting frame defined by the Q"""
      QV = Qmake_a_point(V)
      QV = self.conjugate() * QV * self
      return Vector(QV.q1, QV.q2, QV.q3)
   def set_values(self, V, angle):
      """Sets quaterion values using a direction vector and a rotation of the coordinate frame about it."""
      S = math.sin(-angle/2.)
      self.q1 = V.x * S
      self.q2 = V.y * S
      self.q3 = V.z * S
      self.q4 = math.cos(angle/2.)

   def set_as_QX(self, angle):
      """Sets quaterion in place like QX function"""
      self.q1 = math.sin(-angle/2.)
      self.q2 = 0.
      self.q3 = 0.
      self.q4 = math.cos(angle/2.)

   def set_as_QY(self, angle):
      """Sets quaterion in place like QY function"""
      self.q1 = 0.
      self.q2 = math.sin(-angle/2.)
      self.q3 = 0.
      self.q4 = math.cos(angle/2.)

   def set_as_QZ(self, angle):
      """Sets quaterion in place like QZ function"""
      self.q1 = 0.
      self.q2 = 0.
      self.q3 = math.sin(-angle/2.)
      self.q4 = math.cos(angle/2.)

   def set_as_mult(self, QQ1, QQ2):
      """Sets self as QQ1*QQ2 in place for quaternion multiplication """
      self.q1 = QQ2.q1 * QQ1.q4 + QQ1.q1 * QQ2.q4 + (QQ1.q2 * QQ2.q3 - QQ1.q3 * QQ2.q2)
      self.q2 = QQ2.q2 * QQ1.q4 + QQ1.q2 * QQ2.q4 + (QQ1.q3 * QQ2.q1 - QQ1.q1 * QQ2.q3)
      self.q3 = QQ2.q3 * QQ1.q4 + QQ1.q3 * QQ2.q4 + (QQ1.q1 * QQ2.q2 - QQ1.q2 * QQ2.q1)
      self.q4 = QQ1.q4 * QQ2.q4 - (QQ1.q1 * QQ2.q1 + QQ1.q2 * QQ2.q2 + QQ1.q3 * QQ2.q3)

   def set_as_point(self, V):
      self.q1 = V.x
      self.q2 = V.y
      self.q3 = V.z
      self.q4 = 0.

   def set_equal(self, Q):
      """Assigns values from other Q to this one. """
      self.q1 = Q.q1
      self.q2 = Q.q2
      self.q3 = Q.q3
      self.q4 = Q.q4

   def set_as_conjugate(self):
      """Assigns conjugate values in place. """
      self.q1 *= -1.
      self.q2 *= -1.
      self.q3 *= -1.


class NumericList (list):
	"""List class that supports multiplication.  Only valid for numbers."""
	
	def __mul__ (L1, L2):
		"""Take the dot product of two numeric lists.
		Not using Vector for this because it is limited to three dimensions.
		Lists must have the same number of elements."""
		
		return(sum(map(lambda x, y: x*y, L1, L2)))


class Matrix (list):	
    """Class to encapsulate matrix data and methods.

   	A matrix is simply a list of lists that correspond to rows of the matrix.
	This is just intended to handle simple multiplication and vector rotations.
	For anything more advanced or computationally intensive, Python library routines
	should be used."""

    def __init__(self, rows):
        """Constructor for a matrix.

        This accepts a list of rows.
        It is assumed the rows are all of the same length."""

        for row in rows:
          self.append(NumericList(row))  #copy list	
    	
    def __str__(self):
        """Returns a string representation of the matrix."""

        return_str = 'Matrix:'

        for row_index in range(len(self)):
            row_str = 'Row %d: ' %(row_index + 1)
            row = self[row_index]

            for col_index in range(len(row)):
                row_str = row_str + '%6.3f  ' % (row[col_index])

            return_str = return_str + '\n' + row_str

        return(return_str)

    def element(self, row_index, col_index):
        """Returns an element of the matrix indexed by row and column.

        Indices begin with 0."""
		
        return ((self[row_index])[col_index])
	
    def row(self, row_index):
        """Returns a specified row of the matrix."""
	
        return(self[row_index])

    def column(self, col_index):
        """Returns a specified column of the matrix as a numeric list."""

        return(NumericList([row[col_index] for row in self]))

    def num_rows(self):
        """Returns the number of rows in the matrix."""

        return(len(self))

    def num_cols(self):
        """Returns the number of columns in the matrix."""

        return (len(self[0]))  #assumes all rows of equal length

    def get_cols (self):
        """Returns list of all columns in a matrix."""

        return ([self.column(col_index) for col_index in range(0, self.num_cols())])
	
    def __mul__(m1, m2):
        """Multiplies two Matrix objects and returns the resulting matrix.

        Number of rows in m1 must equal the number of columns in m2."""

        result_rows = []

        #Iterate over the rows in m1.  The first column of row i is formed by
        #multiplying the ith row of m1 by the first column of m2.  The second
        #column is formed by muliplying the ith row of m1 by the second column
        #of m2, etc.
        for row in m1:
        	new_row = []
        	
        	for col in m2.get_cols():
        		new_row.append(row * col)
        		
        	result_rows.append(new_row)

        return (Matrix(result_rows))

	
class Vector (object):
   "Class to encapsulate vector data and operations."

   def __init__(self, x=0.0, y=0.0, z=0.0):
      """Constructor for a three-dimensional vector.

	  Note that two-dimensional vectors can be constructed by omitting one of
	  the coordinates, which will default to 0."""
		
      self.x = x     #Cartesian x coordinate
      self.y = y     #Cartesian y coordinate
      self.z = z     #Cartesian z coordinate

   def __str__(self):  #Called when used in print statement
   	  """Returns a string representation of the vector."""
   	  return('Vector: x: %.3f, y: %.3f, z: %.3f'\
   	  			% (self.x, self.y, self.z)) 	

   def set_eq(self, x=None, y=None, z=None):
      """Assigns new value to vector.

      Arguments are now optional to permit this to be used with 2D vectors
      or to modify any subset of coordinates."""

      if (x != None):
      	self.x = x
      if (y != None):
        self.y = y
      if (z != None):
          self.z = z

   def length(self):
      """Returns magnitude of the vector """
      return(math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z))

   def normalize(self):
      """Returns copy of the normalized vector """
      mag = self.length()
      return (Vector(self.x/mag, self.y/mag, self.z/mag))

   def __mul__(self, rs):
      """Implements Vector * scalar.  Can then use '*' syntax in multiplying a vector by a scalar rs. """
      x = self.x * rs
      y = self.y * rs
      z = self.z * rs
      return (Vector(x, y, z))

   def __rmul__(self, ls):
      """Implements float * Vector """
      x = self.x * ls
      y = self.y * ls
      z = self.z * ls
      return (Vector(x, y, z))

   def __add__(self, rs):
      """Implements Vector + Vector """
      x = self.x + rs.x
      y = self.y + rs.y
      z = self.z + rs.z
      return (Vector(x, y, z))

   def __sub__(self, rs):
      """Implements Vector - Vector """
      x = self.x - rs.x
      y = self.y - rs.y
      z = self.z - rs.z
      return (Vector(x, y, z))

   def __truediv__(self, rs):
      """Implements Vector / float """
      x = self.x / rs
      y = self.y / rs
      z = self.z / rs
      return (Vector(x, y, z))

   def __imul__(self, rs):
      """Implements Vector *= float """
      self.x *= rs
      self.y *= rs
      self.z *= rs
      return (self)

   def __iadd__(self, rs):
      """Implements Vector += vector """
      self.x += rs.x
      self.y += rs.y
      self.z += rs.z
      return (self)

   def __isub__(self, rs):
      """Implements Vector -= vector """
      self.x -= rs.x
      self.y -= rs.y
      self.z -= rs.z
      return (self)

   def __idiv__(self, rs):
      """Implements Vector /= float """
      self.x /= rs
      self.y /= rs
      self.z /= rs
      return (self)

   def create_matrix(self):
       """Converts a Vector into a single-column matrix."""

       column = [self.x, self.y, self.z]
       return(Matrix([[element] for element in column]))  #singleton list

   #Recommend deletion -- better to use a single interface that takes two vectors.
   def dot(self, V2):
      """returns dot product between two vectors """
      return(self.x * V2.x + self.y * V2.y + self.z * V2.z)

   #Recommend deletion in favor of non-method version.
   def cross(self, V2):
       """returns cross product of two vectors """
       x = self.y*V2.z - V1.z*V2.y
       y = self.z*V2.x - V1.x*V2.z
       z = self.x*V2.y - V1.y*V2.x
       return Vector(x, y, z)

   #Replace by separation - RLH
   def angle(self, V2):
      """returns angle between the two vectors in degrees """
      R1 = self.length()
      R2 = V2.length()
      adot = dot(self, V2)
      adot = adot / R1 / R2
      adot = min(1., adot)
      adot = max(-1., adot)
      return acosd(adot)

   # RLH: What do these add?  We're creating methods just to access individual attributes.
   def rx(self):  return self.x
   def ry(self):  return self.y
   def rz(self):  return self.z

   # RLH: Suggest deletion in favor of __str__, which has the advantage that it is called on print.
   def display(self):
      return "[%f, %f, %f]" % (self.x, self.y, self.z)

   # RLH: Not necessary if CelestialVector is used.
   def set_xyz(self, ra, dec):
      """Creates a unit vector from spherical coordinates """
      self.x = cosd(dec) *cosd(ra)
      self.y = cosd(dec) *sind(ra)
      self.z = sind(dec)


class CelestialVector (Vector):
    "Class to encapsulate a unit vector on the celestial sphere."
	
    def __init__(self, ra=0.0, dec=0.0, frame='eq', degrees=True):
        """Constructor for a celestial vector.
		
	    There are two spherical coordinates, a longitudinal coordinate (called
	    right ascension), and a latitudinal coordinate (called declination).
	    The RA is defined as the counterclockwise angle from a reference direction
	    on the equatorial plane; it ranges from 0-360 degrees.  The DEC is the angle
	    between the vector and the equatorial plane; it ranges from -90 to 90 degrees.
	    Angles are specified in degrees but represented internally as radians.
		
	    The frame attribute indicates the coordinate frame of the vector, which may be
	    'eq' (equatorial, default), 'ec' (ecliptic), or 'gal' (galactic).  In equatorial
	    coordinates, the equatorial plane is the celestial equator (extension of the Earth's
	    equator) and the reference axis is the vernal equinox.  In ecliptic coordiantes,
	    the equatorial plane is the ecliptic (the Earth's orbital plane) and the reference
	    axis is usually defined relative to the Sun.  In galactic coordinates, the equatorial
	    plane is the plane of the Galaxy.
	
	    The degrees attribute should be True if the RA, DEC inputs are in degrees.
	    Otherwise radians is assumed.
		
	    The coordinates "ra" and "dec" may be used in all three systems.  Other names for
	    coordinates in different frames may be defined for clarity.
		
	    A CelestialVector is also an ordinary unit vector, with Cartesian coordinates defined
	    relative to the equatorial plane."""
	
        if (degrees):
          ra = math2.D2R * ra
          dec = math2.D2R * dec

        self.ra = ra
        self.dec = dec
        self.frame = frame

        #Initialize standard vector with translated Cartesian coordinates
        Vector.__init__(self, x=cos(ra)*cos(dec), y=sin(ra)*cos(dec), z=sin(dec))
		
    def __str__(self, verbose=True):
      """Returns a string representation of the vector.  Displays angles in degrees."""
      celest_info = 'CelestialVector: RA: %.3fD, DEC: %.3fD, frame: %s'\
      % (math2.R2D*self.ra, math2.R2D*self.dec, self.frame)

      if (verbose):
          celest_info = celest_info + '\n' + super(CelestialVector, self).__str__()
      return(celest_info)
	
    def set_eq(self, ra, dec, degrees=False):
      """Modifies a celestial vector with a new RA and DEC.

      degrees = True if units are degrees.  Default is radians."""

      if (degrees):
          ra = math2.D2R * ra
          dec = math2.D2R * dec

      self.ra = ra
      self.dec = dec

      #Update Cartesian coordinates as well.
      super(CelestialVector, self).set_eq(cos(ra)*cos(dec), sin(ra)*cos(dec), sin(dec))
		
    def update_cartesian(self, x=None, y=None, z=None):
      """Modifies a celestial vector by specifying new Cartesian coordinates.

         Any subset of the Cartesian coordinates may be specifed."""

      if (x != None):
         self.x = x
      if (y != None):
         self.y = y
      if (z != None):
       	 self.z = z
      	
      self.ra = atan2(self.y, self.x)  #RA is arctan of y/x
      if (self.ra < 0):				 #Make sure RA is positive
         self.ra += 2*pi
      	
      self.dec = math2.asin2(self.z)  #DEC is arcsin of z
		
    def rotate_about_axis (self, angle, axis):
        """This rotates a vector about an axis by the specified angle
        by using a rotation matrix.
        A new vector is returned.

        Axis must be 'x', 'y', or 'z'.
        The x-rotation rotates the y-axis toward the z-axis.
        The y-rotation rotates the z-axis toward the x-axis.
        The z-rotation rotates the x-axis toward the y-axis."""

        if (axis == 'x'):
            rot_matrix = Matrix([[1, 0, 0], [0, cos(angle), -sin(angle)], \
            [0, sin(angle), cos(angle)]])

        elif (axis == 'y'):
        	rot_matrix = Matrix([[cos(angle), 0, sin(angle)], [0, 1, 0], \
        	[-sin(angle), 0, cos(angle)]])
        	
        elif (axis == 'z'):
        	rot_matrix = Matrix([[cos(angle), -sin(angle), 0], \
        	[sin(angle), cos(angle), 0], [0, 0, 1]])
        	
        else:
            print ('Error')
            return

        new_matrix = rot_matrix * self.create_matrix()
        new_vector = new_matrix.column(0)
        result = CelestialVector()  #initialize with Cartesian coordiantes
        result.update_cartesian(x=new_vector[0], y=new_vector[1], z=new_vector[2])
        return(result)

    def rotate_about_eigenaxis(self, angle, eigenaxis):
        """rotates a vector about arbitrary eigenaxis.

        eigenaxis = Vector object (axis about which to rotate).
        angle = angle to rotate by in radians.
        Rotation is counterclockwise looking outward from origin along eigenaxis.
        Function uses rotation matrix from Rodrigues formula.

        Note: This function is more general than rotate_about_axis above and
        could be used in its place.  However, rotate_about_axis is faster and
        clearer when the rotation axis is one of the Cartesian axes."""

        cos_ang = cos(angle)    #Used repeatedly below
        sin_ang = sin(angle)

        #Fill out the Rodrigues rotation matrix
        R11 = cos_ang + eigenaxis.x**2 * (1 - cos_ang)
        R12 = eigenaxis.x * eigenaxis.y * (1 - cos_ang) - eigenaxis.z * sin_ang
        R13 = eigenaxis.x * eigenaxis.z * (1 - cos_ang) + eigenaxis.y * sin_ang
        R21 = eigenaxis.x * eigenaxis.y * (1 - cos_ang) + eigenaxis.z * sin_ang
        R22 = cos_ang + eigenaxis.y**2 * (1 - cos_ang)
        R23 = eigenaxis.y * eigenaxis.z * (1 - cos_ang) - eigenaxis.x * sin_ang
        R31 = eigenaxis.x * eigenaxis.z * (1 - cos_ang) - eigenaxis.y * sin_ang
        R32 = eigenaxis.y * eigenaxis.z * (1 - cos_ang) + eigenaxis.x * sin_ang
        R33 = cos_ang + eigenaxis.z**2 * (1 - cos_ang)

        rot_matrix = Matrix([[R11, R12, R13], [R21, R22, R23], [R31, R32, R33]])
        new_matrix = rot_matrix * self.create_matrix()
        new_vector = new_matrix.column(0)
        result = CelestialVector()  #initialize with Cartesian coordinates
        result.update_cartesian(x=new_vector[0], y=new_vector[1], z=new_vector[2])
        return(result)

    def rotate_using_quaternion(self, angle, eigenaxis):
        """Rotates a vector about arbitrary eigenaxis using quaternion.

        This is an alternative formulation for rotate_about_eigenaxis.
        Interface is the same as rotate_about_eigenaxis."""

        q = Quaternion(eigenaxis, 0.0)

        #Need to negate here because set_values performs a negative rotation
        q.set_values(eigenaxis, -angle)   #quaternion now represents the rotation
        return(make_celestial_vector(q.cnvrt(self)))

    def transform_frame(self, new_frame):
        """Transforms coordinates between celestial and ecliptic frames

        and returns result as a new CelestialVector.
        If new coordinate frame is the same as the old, a copy of the vector
        is returned."""

        result = None

        #Equatorial to ecliptic: rotate z-axis toward y-axis.
        if ((new_frame == 'ec') and (self.frame == 'eq')):
            result = self.rotate_about_axis(-math2.OBLIQUITY, 'x')

        #Ecliptic to equatorial: rotate y-axis toward z-axis.
        elif ((new_frame == 'eq') and (self.frame == 'ec')):
        	result = self.rotate_about_axis(math2.OBLIQUITY, 'x')
        	
        elif ((new_frame == 'gal') and (self.frame == 'eq')):
            #Use formula from Wayne Kinzel's book, adjusted for J2000 coordinates.
            b = math2.asin2(cos(self.dec) * cos(NGP.longitude) * cos(self.ra - NGP.latitude)\
            + sin(self.dec) * sin(NGP.longitude))

            l = atan2(sin(self.dec) - sin(b) * sin (NGP.longitude), \
            cos (self.dec) * sin(self.ra - NGP.latitude) * cos(NGP.longitude)) + NGP.anode

            result = CelestialVector(l, b, degrees=False)

        elif ((new_frame == 'eq') and (self.frame == 'gal')):
            l = self.ra   #use l, b notation here for clarity
            b = self.dec

            dec = math2.asin2(cos(b) * cos(NGP.longitude) * sin(l - NGP.anode) + sin(b) * sin(NGP.longitude))

            ra = atan2(cos(b) * cos(l - NGP.anode), \
            sin(b) * cos(NGP.longitude) - cos(b) * sin(NGP.longitude) * sin(l - NGP.anode)) + NGP.latitude

            result = CelestialVector(ra, dec, degrees=False)
                     	
        elif (((new_frame == 'gal') and (self.frame == 'ec')) or ((new_frame == 'ec') and (self.frame == 'gal'))):
        	print("Error: Direct conversion between ecliptic and galactic coordinates not supported yet")
        	
        elif (new_frame != self.frame):
        	print("Error: unrecognized coordinate frame.")

        #If there was an error, return a copy of the initial vector.
        if (result is None):
            result = CelestialVector(self.ra, self.dec, self.frame, False)

        else:
            result.frame = new_frame  #record new frame

        return (result)

    def rotate_by_posang (self, pa):
        """Returns the vector that results from rotating the self vector
        counterclockwise from the North projection onto the plane
        orthogonal to that vector by the specified position angle
        (in radians). See "V3-axis Position Angle", John Isaacs, May 2003 for
        further discussion."""
        x_coord = -cos(self.ra) * sin(self.dec) * cos(pa) - sin(self.ra) * sin(pa)
        y_coord = -sin(self.ra) * sin(self.dec) * cos(pa) + cos(self.ra) * sin(pa)
        z_coord = cos(self.dec) * cos(pa)
        result = CelestialVector()
        result.update_cartesian(x_coord, y_coord, z_coord)
        return(result)

    def position_angle (self, v):
        """Returns the position angle of v at the self vector, in radians.

        v is an arbitrary vector that should be a CelestialVector object.
        The position angle is the angle between the North vector on the
        plane orthogonal to the self vector and the projection of v onto
        that plane, defined counterclockwise.
        See "V3-axis Position Angle", John Isaacs, May 2003 for
        further discussion."""
        y_coord = cos(v.dec) * sin(v.ra - self.ra)
        x_coord = sin(v.dec) * cos(self.dec) - cos(v.dec) * sin(self.dec) * cos (v.ra - self.ra)
        pa = atan2(y_coord, x_coord)

        if (pa < 0):
            pa += (2 * pi)  #PA has range 0-360 degrees

        return(pa)


class Attitude (CelestialVector):
    "Defines an Observatory attitude by adding a position angle."""

    def __init__(self, ra=0.0, dec=0.0, pa=0.0, frame='eq', degrees=True):
        """Constructor for an Attitude.

        pa = position_angle in degrees(default) or radians if degrees=False is specified.
        Other arguments are the same as with CelestialVector."""

        super(Attitude, self).__init__(ra=ra, dec=dec, frame=frame, degrees=degrees)

        if (degrees):   #convert into radians
            pa = math2.D2R * pa

        self.pa = pa

    def __str__ (self, verbose=True):
        """Returns a string representation of the attitude.

        verbose (optional) = flag indicating whether detailed Vector information should be included."""

        att_info = 'Attitude: PA: %.3fD' %(math2.R2D * self.pa)
        att_info = att_info + '\n' + super(Attitude, self).__str__(verbose)
        return(att_info)

#Functions that operate on vectors but are not methods.
def dot(v1, v2):
   """returns dot product between two vectors, non class member """

   return(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

def cross(v1, v2):
   """returns cross product between two vectors, non class member """
   x = v1.y*v2.z - v1.z*v2.y
   y = v1.z*v2.x - v1.x*v2.z
   z = v1.x*v2.y - v1.y*v2.x
   return Vector(x, y, z)

def separation(v1, v2, norm=False):
    """Returns angle between two unit vectors in radians.
   	
    The angle between two normalized vectors is the arc-cosine of the dot product.
    Unless the norm attribute is set to True, it is assumed the vectors are
    already normalized (for performance)."""

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
   	
def ra_delta (v1, v2):
    """Returns difference in right ascension between two CelestialVectors."""

    delta_ra = v1.ra - v2.ra

    # Check for zero crossings.  If the difference exceeds 180 degrees, adjust by
    # 360 in opposite direction.

    if (delta_ra < -pi):
        delta_ra = delta_ra + 2*pi
    elif (delta_ra > pi):
        delta_ra = delta_ra - 2*pi

    return(delta_ra)
   	
def ra_separation(v1, v2):
    """Returns separation in right ascension between two CelestialVectors.
    This is accurate only if the difference in declination is small.

    |sep| = DELTA-RA cos DEC."""

    delta_ra = ra_delta(v1, v2)
    dec = math2.avg2(v1.dec, v2.dec)  #use average of the two declinations.
    return(delta_ra * cos(dec))

def dec_separation(v1, v2):
    """Returns difference in declination between two CelestialVectors."""

    return(v1.dec - v2.dec)    #simply take the difference in declination
   	
def make_celestial_vector(v):
	"""Takes a Vector object and creates an equivalent CelestialVector.
	
	Input vector v must be a unit vector."""
	
	result = CelestialVector()
	result.update_cartesian(v.x, v.y, v.z)
	return(result)
	
def projection (v, axis):
	"""Returns projection of vector v on plane normal to axis.
	
	First take cross-product of v and the axis and normalize it.
	Then cross the axis with the result and return a CelestialVector.
	See http://www.euclideanspace.com/maths/geometry/elements/plane/
    lineOnPlane/index.htm."""
	
	return(make_celestial_vector(cross(axis, (cross(v, axis)).normalize())))

def pos_V_to_ra_dec(V):
   """Returns tuple of spherical angles from unit direction Vector """
   ra = atan2d(V.y, V.x)
   V.z = min(1., V.z)
   V.z = max(-1., V.z)
   dec = asind(V.z)
   if ra < 0.:
      ra += 360.
   return(ra, dec)

# RLH: Recommend replacement by separation.
def angle(V1, V2):
    """returns angle between two vectors in degrees, non class member """
    R1 = V1.length()
    R2 = V2.length()
    adot = dot(V1, V2)
    adot = adot / R1 / R2
    adot = min(1., adot)
    adot = max(-1., adot)
    return acosd(adot)


def vel_ab(U, Vel):
    """ Takes a unit vector and a velocity vector(km/s) and returns a unit
    vector modidifed by the velocity abberation."""
    c = 2.9979e5  # speed of light in km/s
    Beta = Vel * (1./c)
    rgamma = math.sqrt(1.-dot(Beta, Beta))  # This is 1/gamma
    ubeta = dot(U, Beta)
    b = (1./(1. + ubeta))
    return (U*rgamma + Beta * (1. + (1.-rgamma)*ubeta/dot(Beta, Beta))) * b
