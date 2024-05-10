"""
Created on Fri May  3 12:58:19 2024
Hyperbolic Tiling generator
@author: Benjamin Takacs

This program was made alongside my admission thesis ''Hyperbolic Reflection Groups''.
It generates a tiling of the poincarée disk with n-gons where at each vertex m more n-gons meet and saves it as a png.
"""

import numpy as np
import math
from cmath import sqrt
import matplotlib.pyplot as plt

### Values
""" These are the main values for this program they are as follows
vertices: The amount of vertices each polygon has
pol: The amount of Polygons that share the same vertex
depth: This is how far away from the first polygon we reflect.
"""
print("Please input the amount of verteces of each polygon (p), how many polygons meet at each vertex (q) and how far the iteration should run.")
print("Keep in mind, that p and q have to satisfy 1/p+1/q<1/2 and p,q>2")
readytogo=1
while readytogo:
    vertices = int(input("Vertex per polygon: "))
    pol = int(input("Polygons per vertex: "))
    depth = int(input("Depth of iteration: "))
    if 1/vertices + 1/pol >= 1/2:
        print("p and q don't meet the requirements, please input correct values")
    else:
        readytogo=0
### Classes ###
"""
Here we define the classes.

"""
"""
The class hypPoint represents a vector in the Lorentzian 3-space
Here the coordinates are saved and the vector can tell if it is in H^2 or is 
space-like. We can also give out the Lorentzian norm of the vector.
"""
class hypPoint:                     # In this we define a point in R^3_L
    def __init__(self, x=1,y=0,z=0):
        self.x=x
        self.y=y
        self.z=z
    def isH2(self):                 # To Check if the point is in H^2
        if abs(-self.x**2+self.y**2+self.z**2+1)<10**(-5) and self.x>0:
            return 1
        else: return 0
    def isSpaceLike(self):
        if -self.x**2+self.y**2+self.z**2+1>0:
            return 1
        else: return 0
        
    def norm(self):                 # This gives back the Lorentzian norm
        return sqrt(-self.x**2+self.y**2+self.z**2)
    def __str__(self):              # Gives back the coordinates of the point
        return f"({self.x}, {self.y}, {self.z})"

"""
We use the class pPoints to represent pints in R^2. These are necessary to plot
our tiling.
These points can tell if they are in the conformal ball and their norm
"""
class pPoint:                       # Here we define a Point in R^2
    def __init__(self, x=0,y=0):
        self.x=x
        self.y=y
    def isInSphere(self):             # To see if the point lies in B^2
        if math.sqrt(self.x**2+self.y**2)<1:
            return 1
        else: return 0
    def norm(self):
        return math.sqrt(self.x**2+self.y**2)
    def __str__(self):              # Gives back the coordinates of the point
        return f"({self.x}, {self.y})"

"""
The class geodesic saves the necessary details in order to represent a unique
geodesic. It saves the start and end point, as well as their distance.
Also the orthonormal vector which is necessary to calculate points in the 
geodesic is saved in order.
The orthogonal vector is necessary in order to reflect on the side of a given
geodesic line.
"""
class geodesic:                     # Defines a geodesic of a certain length using orthonormal vectors
    def __init__(self, start=hypPoint(1,0,0), end=hypPoint(1,0,0)):
        self.start=start
        self.end=end
        self.L = hypMetric(start, end)
        variable = hypMult(-np.cosh(self.L), start)
        variable = hypAddition(end, variable)
        variable = hypMult(1/np.sinh(self.L), variable)
        self.orth=variable
        self.orthogonal=orthogonalize(start, end)
    def isOrthonormal(self):         # Checks if the given vectors are orthonormal
        val1 = self.start.norm()**2
        val2 = self.orth.norm()**2
        if abs(val1.real+1)<10**(-5) and abs(val2.real-1)<10**(-5) and abs(lorentzSkalar(self.start, self.orth))<10**(-5):
            return 1
        else: return 0
    def __str__(self):
        return f"({self.start}, {self.end}, {self.orth})"
            
"""
The class polygon consists of the vertices as well as the geodesics between a 
vertex and its neighbours. It also memorizes if it was already used to generate
the neighbour polygons.
"""
class polygon:
    def __init__(self, vertices:[]):
        self.vertices=[]
        for i in vertices:
            if i.isH2():
                self.vertices.append(i)
            else: 
                print("Point " + i + "is not in H2.")
        self.Geodesics = []
        if len(self.vertices)>1:
            for i in range(len(self.vertices)-1):
                self.Geodesics.append(geodesic(self.vertices[i],self.vertices[i+1]))
            self.Geodesics.append(geodesic(self.vertices[-1], self.vertices[0]))
        self.wasReflectedOn=0
    def plotPolygon(self):
        for geod in self.Geodesics:
            plotGeodesic(geod)


""" 
Here we define the functions used to generate the tiling.

"""

### Vector functions ###
"""
Here we have functions necessary to calculate with vectors.
"""
"This gives back the lorentz scalar of two vectors"
def lorentzSkalar(vec1:hypPoint, vec2:hypPoint)-> float: 
    return -vec1.x*vec2.x+vec1.y*vec2.y+vec1.z*vec2.z

"This adds two vectors together"
def hypAddition(vec1:hypPoint, vec2:hypPoint)-> hypPoint:
    return hypPoint(vec1.x+vec2.x,vec1.y+vec2.y,vec1.z+vec2.z)

"This multiplies a vector by some real number"
def hypMult(lambd:float, vec:hypPoint)->hypPoint:       
    return hypPoint(lambd*vec.x, lambd*vec.y, lambd*vec.z)

"This gives back the metric d_H"
def hypMetric(vec1:hypPoint, vec2:hypPoint)-> float:    
    if vec1.isH2() and vec2.isH2():
        return np.arccosh(-lorentzSkalar(vec1, vec2))
    else: return 0
    
"""This gives back a vector orthogonal to two given vectors. It uses the Lorentz
cross product"""
def orthogonalize(vec1:hypPoint, vec2:hypPoint)-> hypPoint: # Generates a vector orthogonal to two given vectors
    x = vec1.z*vec2.y-vec1.y*vec2.z
    y = vec1.z*vec2.x-vec1.x*vec2.z
    z = vec1.x*vec2.y-vec1.y*vec2.x
    return hypPoint(x,y,z)

    
### Generate and compare Polygons ###
"""
Here we have functions that generate and compare Polygons
"""
"""We set the first vertex on a point (x,0). This gives us back the x-koordinate.
Note that we use the conformal ball here, because of its properties in regards
to angles."""
def generateFirstVertice(pol:int, vertices:int)->pPoint:
    a = math.cos(math.pi*(vertices+pol)/vertices/pol) 
    b = math.cos(math.pi*(pol-vertices)/vertices/pol)
    return pPoint(math.sqrt(a / b),0)
    
""""Using the first vertex know the radius of the center polygon. By rotating 
the first vertex we can genereate the others."""
def generateInitialPolygon(pol:int, vertices:int)->[]:
    points = []
    angle = 2*math.pi/vertices
    verticeOne = generateFirstVertice(pol, vertices)
    points.append(projectToH(verticeOne))
    for i in range(1,vertices):
        ang=angle*i
        x = math.cos(ang)*verticeOne.x
        y = math.sin(ang)*verticeOne.x
        points.append(projectToH(pPoint(x,y)))
    return points

"""This code is not optimized. But in order to reduce the times we can compare
two polygons if they are equal."""
def isNew(pol1:polygon, pol2:polygon):
    matches = 0
    for i in range(len(pol1.vertices)):
        for j in range(i-1,len(pol2.vertices)):
            xDiff = pol1.vertices[i].x-pol2.vertices[j].x
            yDiff = pol1.vertices[i].y-pol2.vertices[j].y
            zDiff = pol1.vertices[i].z-pol2.vertices[j].z
            if abs(xDiff)<10**(-2) and abs(yDiff)<10**(-2) and abs(zDiff)<10**(-2):
                matches = matches + 1 
    if matches == len(pol1.vertices):
        return 0
    else:
        return 1
### Reflections ###
"""
Here we have functions that generate reflections. These are used to generate the
tiling.
"""
"This reflects a single vector on the hyperplane P_v"
def reflect(x:hypPoint, v:hypPoint)->hypPoint:
    mod = -2*lorentzSkalar(x, v)/lorentzSkalar(v, v)
    var = hypMult(mod, v)
    val = hypAddition(x, var)
    if not v.isSpaceLike():
       print("Warnung") 
    return val

"Using the reflect function we can reflect every point of a polygon on one hyperplane."
def reflectPolygon(pol:polygon, v:hypPoint)->polygon:
    if (v.norm()**2).real>0:
        points = pol.vertices
        nPoints =[]
        for i in points:
            nPoints.append(reflect(i, v))
        return polygon(nPoints)

"Here we reflect a polygon on every of its sides."    
def reflectAllSides(pol:polygon)->[]:
    Polygons = []
    for geod in pol.Geodesics:
        Polygons.append(reflectPolygon(pol, geod.orthogonal))
    pol.wasReflectedOn=1
    return Polygons

"""This is the heart of this program. It generates the first polygon and uses
reflections to generate the adjacent polygons. By checking if a polygon was reflected
already and if a reflected polygon is new, we can reduce the processing time.
It returns the array of all the polygons that will be plotted."""
def iterate(vertices:int, polPerVert:int, depth:int)->[]:
    Polygons = [polygon(generateInitialPolygon(pol, vertices))]
    for i in range(depth):
        pols = []
        for p in Polygons:
            if p.wasReflectedOn == 0:  
                newPols = reflectAllSides(p)
                for i in newPols:
                    pols.append(i)
        for newPol in pols:
            for oldPol in Polygons:
                if not isNew(newPol,oldPol):
                    try:
                        pols.remove(newPol)
                    except:pass
        for i in pols:
            Polygons.append(i)
    return Polygons
    
### Projections and Visualization ###
"""
With these functions we use the stereographic projection and plot a geodesic.
"""
"This returns the projection of a vector in H^2 to B^2."
def projectToB(vec:hypPoint)-> pPoint:  # Projects a Vector in H^2 on B^2
    if vec.isH2():
        try:
            output = pPoint(vec.y/(1+vec.x), vec.z/(1+vec.x))
            return output
        except ZeroDivisionError: print("Nicht durch 0 Teilen!")
    else: print("Kein Element von H^n!")
    
"This returns the projection of a vector in B^2 to H^2."
def projectToH(vec:pPoint())-> hypPoint:  # Projects a Vector in H^2 on B^2
    if vec.isInSphere():
        try:
            mod = 1/(1-vec.norm()**2)
            return hypPoint((1+vec.norm()**2)*mod,2*vec.x*mod,2*vec.y*mod)
        except ZeroDivisionError: print("Nicht durch 0 Teilen!")
    else: print("Kein Element von B^2!")
    
"This returns the value of a geodesic on a given input t as a pPoint"
def geodesicPoint(gam:geodesic,t:float)-> pPoint:   # Returns the value of some geodesic gam at a point t
    if t<=gam.L and t>=0 and gam.isOrthonormal():
        nVec1 = hypMult(np.cosh(t), gam.start)
        nVec2 = hypMult(np.sinh(t), gam.orth)
        return projectToB(hypAddition(nVec1, nVec2))
    
"""Using the geodesicPoint funtion we can plot the geodesic by dividing it into 
100 points"""
def plotGeodesic(geod:geodesic): # Plots a geodesic Line in B2
    xs = []
    ys = []
    for i in range(100):
        t = i*geod.L/100
        point = geodesicPoint(geod, t)
        xs.append(point.x)
        ys.append(point.y)
    plt.plot(xs,ys, color='black')
    
### Main ###
"Here we generate the polygons and plot them."
Pols = iterate(vertices,pol,depth)
for i in Pols:
    i.plotPolygon()    

"Here we generate a unit circle as well as the coordinate axis"
angles = np.linspace(0 * math.pi, 2 * math.pi, 100) 
xs = np.cos(angles)
ys = np.sin(angles)

plt.plot(xs, ys, color = 'black')
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)
plt.gca().set_aspect('equal')
plt.draw()
name = "Tiling-"+str(vertices)+"-"+str(pol)+"-"+str(depth)+".png"
plt.savefig(name, format='png')