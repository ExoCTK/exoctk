#! /usr/bin/env python
#quaternion module

"""Version 4 September 9, 2010 WMK
Flipped sign of the angle in the QX,QY,QZ,QJX,QJY,QJZ, set_values, set_as_QX,...
functions to be consistent with the corrected multiplication. Also updated
the doc strings."""

"""Version 3 September 8, 2010 RLH
Backed out change to cnvrt in version 2."""

"""Version 2 September 3, 2010 RLH
Fixed sign error in quaternion multiplication functions.
Added quaternion __str__ method.
Modified cnvrt method to return a CelestialVector."""


#4/9/2010 WMK
#Redefined the __init__ inputs
#Changed from a Vector internal representation to 3 scalers
#Fixed an error in cvt_att_Q_to_angles, was assuming an att2inertial Quaternion!
#Streamlined some of the functions
#  Version 1.0 August 3, 2010
#  Got rid of degrees trig functions.

from math import *
from ExoCTK.tor.tor2.math_extensionsx import *
# import ExoCTK.tor.contam_tool.rotationsx as ROT

D2R = pi/180.
R2D = 180. / pi
PI2 = 2. * pi
unit_limit = lambda x: min(max(-1.,x),1.)




# This object is from rotationsx.py. I had to copy it here because DL is doing cyclic imports
class Vector (object):
   "Class to encapsulate vector data and operations."
     
   def __init__(self,x=0.0,y=0.0,z=0.0):
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
      return(sqrt(self.x * self.x + self.y * self.y + self.z * self.z))
       
   def normalize(self):
      """Returns copy of the normalized vector """ 
      mag = self.length()
      return (Vector(self.x/mag,self.y/mag,self.z/mag))

   def __mul__(self,rs):
      """Implements Vector * scalar.  Can then use '*' syntax in multiplying a vector by a scalar rs. """ 
      x = self.x * rs
      y = self.y * rs
      z = self.z * rs
      return (Vector(x,y,z))
      
   def __rmul__(self,ls):
      """Implements float * Vector """ 
      x = self.x * ls
      y = self.y * ls
      z = self.z * ls
      return (Vector(x,y,z))
      
   def __add__(self,rs):
      """Implements Vector + Vector """ 
      x = self.x + rs.x
      y = self.y + rs.y
      z = self.z + rs.z
      return (Vector(x,y,z))
      
   def __sub__(self,rs):
      """Implements Vector - Vector """ 
      x = self.x - rs.x
      y = self.y - rs.y
      z = self.z - rs.z
      return (Vector(x,y,z))
      
   def __truediv__(self,rs):
      """Implements Vector / float """ 
      x = self.x / rs
      y = self.y / rs
      z = self.z / rs
      return (Vector(x,y,z))
      
   def __imul__(self,rs):
      """Implements Vector *= float """ 
      self.x *= rs
      self.y *= rs
      self.z *= rs
      return (self)
      
   def __iadd__(self,rs):
      """Implements Vector += vector """ 
      self.x += rs.x
      self.y += rs.y
      self.z += rs.z
      return (self)
      
   def __isub__(self,rs):
      """Implements Vector -= vector """ 
      self.x -= rs.x
      self.y -= rs.y
      self.z -= rs.z
      return (self)
      
   def __idiv__(self,rs):
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
   def dot(self,V2):
      """returns dot product between two vectors """ 
      return(self.x * V2.x + self.y * V2.y + self.z * V2.z)

   #Recommend deletion in favor of non-method version.
   def cross(self,V2):
       """returns cross product of two vectors """ 
       x = self.y*V2.z - V1.z*V2.y
       y = self.z*V2.x - V1.x*V2.z
       z = self.x*V2.y - V1.y*V2.x
       return Vector(x,y,z)
      
   #Replace by separation - RLH
   def angle(self,V2):
      """returns angle between the two vectors in degrees """ 
      R1 = self.length()
      R2 = V2.length()
      adot = dot(self,V2)
      adot = adot / R1 / R2
      adot = min(1.,adot)
      adot = max(-1.,adot)
      return acosd(adot)
            
   # RLH: What do these add?  We're creating methods just to access individual attributes.
   def rx(self):  return self.x
   def ry(self):  return self.y
   def rz(self):  return self.z
   
   # RLH: Suggest deletion in favor of __str__, which has the advantage that it is called on print.
   def display(self):
      return "[%f, %f, %f]" % (self.x,self.y,self.z)
      
   # RLH: Not necessary if CelestialVector is used.
   def set_xyz(self,ra,dec):
      """Creates a unit vector from spherical coordinates """ 
      self.x = cosd(dec) *cosd(ra)
      self.y = cosd(dec) *sind(ra)
      self.z = sind(dec)
#End of the vector object




def QX(angle):
   """Creates rotation quaternion about X axis, rotates a vector about this axis"""
   return Quaternion(Vector(sin(angle/2.),0.,0.),cos(angle/2.))

def QY(angle):
   """Creates rotation quaternion about Y axis, rotates a vector about this axis"""
   return Quaternion(Vector(0.,sin(angle/2.),0.),cos(angle/2.))

def QZ(angle):
   """Creates rotation quaternion about Z axis, rotates a vector about this axis"""
   return Quaternion(Vector(0.,0.,sin(angle/2.)),cos(angle/2.))

def Qmake_a_point(V):
   """Creates a pure Q, i.e. defines a pointing not a rotation"""
   return Quaternion(V,0.)

def cvt_pt_Q_to_V(Q):
   """Converts a pure (pointing) Q to a unit position Vector"""
   return Vector(Q.q1,Q.q2,Q.q3)


#The following functions are dependent upon the spacecraft definitions and perhaps should be moved to that module

def Qmake_body2inertial(coord1,coord2,V3pa):
   """Creates a rotation Q, going from the body frame to inertial"""
   return QZ(coord1)*QY(-coord2)*QX(-V3pa)

def Qmake_v2v3_2body(v2,v3):
   """Creates a rotation Q, going from v2 and v3 in the body frame to inertial"""
   return QY(v3)*QZ(-v2)

def Qmake_v2v3_2inertial(coord1,coord2,V3pa,v2,v3):
   """Creates a rotation Q, going from v2 and v3 in the body frame to inertial"""
   return QZ(coord1)*QY(-coord2)*QX(-V3pa)*QY(v3)*QZ(-v2)

def Qmake_aperture2inertial(coord1,coord2,APA,xoff,yoff,s,YapPA,V3ref,V2ref):
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
   coord1 = atan2(r21,r11)
   if coord1 < 0. : coord1 += PI2
   coord2 = asin2(r31)  # use "safe" version of sine
   pa  = atan2(-r32,r33)
   if pa < 0. : pa += PI2
   return (coord1,coord2,pa)

def cvt_v2v3_using_body2inertial_Q_to_c1c2pa_tuple(Q,v2,v3):
   """Given Q and v2,v3 gives pos on sky and V3 PA """
   Vp_body = Vector(0.,0.,0.)
   Vp_body.set_xyz_from_angs(v2,v3)
   Vp_eci_pt = Q.cnvrt(Vp_body)
   coord1  = atan2(Vp_eci_pt.y,Vp_eci_pt.x)
   if coord1 < 0. : coord1 += PI2
   coord2 = asin(unit_limit(Vp_eci_pt.z))

   V3_body = Vector(0.,0.,1.)
   V3_eci_pt = Q.cnvrt(V3_body)
   NP_eci = Vector(0.,0.,1.)
   V_left = cross(NP_eci,Vp_eci_pt)
   if V_left.length()>0.:
      V_left = V_left/V_left.length()
   NP_in_plane = cross(Vp_eci_pt,V_left)
   x = dot(V3_eci_pt,NP_in_plane)
   y = dot(V3_eci_pt,V_left)
   pa  = atan2(y,x)
   if pa < 0. : pa += PI2

   return (coord1,coord2,pa)

def cvt_c1c2_using_body2inertial_Q_to_v2v3pa_tuple(Q,coord1,coord2):
   """Given Q and a position, returns v2,v3,V3PA tuple """
   Vp_eci = Vector(1.,0.,0.)
   Vp_eci.set_xyz_from_angs(coord1,coord2)
   Vp_body_pt = Q.inv_cnvrt(Vp_eci)
   v2 = atan2(Vp_body_pt.y,Vp_body_pt.x)
   v3 = asin(unit_limit(Vp_body_pt.z))
   V3_body = Vector(0.,0.,1.)
   V3_eci_pt = self.cnvrt(V3_body)
   NP_eci = Vector(0.,0.,1.)
   V_left = cross(NP_eci,Vp_eci)
   if V_left.length()>0.:
    V_left = V_left / V_left.length()
   NP_in_plane = cross(Vp_eci,V_left)
   x = dot(V3_eci_pt,NP_in_plane)
   y = dot(V3_eci_pt,V_left)
   pa  = atan2(y,x)
   if pa < 0. : pa += PI2
   return (v2,v3,pa)


############################################################


class Quaternion:
   """This representation is used by Wertz and Markley """
   def __init__(self,V,q4):
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
      return math.sqrt(self.q1*self.q1 + self.q2*self.q2 + self.q3*self.q3 + self.q4*self.q4)
   
   def normalize(self):
      """Returns a copy of the Q normalized """
      scale = self.length()
      return Quaternion(Vector(self.q1/scale,self.q2/scale,self.q3/scale),self.q4/scale)
   
   def conjugate(self):
      """Returns a copy of the conjugated Q """
      return Quaternion(Vector(-self.q1,-self.q2,-self.q3),self.q4)
   
   def __mul__(self,rs):
      """Defines Q*Q for quaternion multiplication """
      Q = Quaternion(Vector(0.,0.,0.),0.)
      #Q.V = rs.V * self.q4 + self.V * rs.q4 + cross(self.V,rs.V)
      Q.q1 = rs.q1 * self.q4 + self.q1 * rs.q4 + (self.q2 * rs.q3 - self.q3 * rs.q2)
      Q.q2 = rs.q2 * self.q4 + self.q2 * rs.q4 + (self.q3 * rs.q1 - self.q1 * rs.q3)
      Q.q3 = rs.q3 * self.q4 + self.q3 * rs.q4 + (self.q1 * rs.q2 - self.q2 * rs.q1)
      Q.q4 = self.q4 * rs.q4 - (self.q1 * rs.q1 + self.q2 * rs.q2 + self.q3 * rs.q3)
      return Q
      
   def cnvrt(self,V):
      """Rotates a vector from the starting frame to the ending frame defined by the Q """
      QV = Qmake_a_point(V)
      QV = self * QV * self.conjugate()
      return Vector(QV.q1,QV.q2,QV.q3)
      
   def inv_cnvrt(self,V):
      """Rotates a vector from the ending frame to the starting frame defined by the Q"""
      QV = Qmake_a_point(V)
      QV = self.conjugate() * QV * self
      return Vector(QV.q1,QV.q2,QV.q3)
   def set_values(self, V, angle):
      """Sets quaterion values using a direction vector and a rotation of the coordinate frame about it."""
      S = sin(-angle/2.)
      self.q1 = V.x * S
      self.q2 = V.y * S
      self.q3 = V.z * S
      self.q4 = cos(angle/2.)
      
   def set_as_QX(self,angle):
      """Sets quaterion in place like QX function"""
      self.q1 = sin(-angle/2.)
      self.q2 = 0.
      self.q3 = 0.
      self.q4 = cos(angle/2.)
      
   def set_as_QY(self,angle):
      """Sets quaterion in place like QY function"""
      self.q1 = 0.
      self.q2 = sin(-angle/2.)
      self.q3 = 0.
      self.q4 = cos(angle/2.)
      
   def set_as_QZ(self,angle):
      """Sets quaterion in place like QZ function"""
      self.q1 = 0.
      self.q2 = 0.
      self.q3 = sin(-angle/2.)
      self.q4 = cos(angle/2.)
      
   def set_as_mult(self,QQ1,QQ2):
      """Sets self as QQ1*QQ2 in place for quaternion multiplication """
      self.q1 = QQ2.q1 * QQ1.q4 + QQ1.q1 * QQ2.q4 + (QQ1.q2 * QQ2.q3 - QQ1.q3 * QQ2.q2)
      self.q2 = QQ2.q2 * QQ1.q4 + QQ1.q2 * QQ2.q4 + (QQ1.q3 * QQ2.q1 - QQ1.q1 * QQ2.q3)
      self.q3 = QQ2.q3 * QQ1.q4 + QQ1.q3 * QQ2.q4 + (QQ1.q1 * QQ2.q2 - QQ1.q2 * QQ2.q1)
      self.q4 = QQ1.q4 * QQ2.q4 - (QQ1.q1 * QQ2.q1 + QQ1.q2 * QQ2.q2 + QQ1.q3 * QQ2.q3)
      
   def set_as_point(self,V):
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

      
   
