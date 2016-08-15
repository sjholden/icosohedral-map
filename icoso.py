# MIT License
# 
# Copyright (c) 2016 Sam Holden
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Create an icosohedral projection of a map.

There is obviously a way to do this with actual math, but this is pure slow brute force.

A gnomonic projection is used for each face of the icosohedron.

The terrible algorithm is:

Create an image containing the triangle mesh. It will 5.5 triangle side lengths wide and 3 triangle heights high.
Iterate over every point in the image.
If the point is inside one of the triangles (we use the simple check if the point is on the same side of the defined by each vertex pair as the remaining vertex):
    Determine the offset of the point within the triangle, the x offset is in units of the side length, the y offset is in units of the triangle height.
    Find the corresponding triangle in an icosohedron.
    Convert the offsets - the x offset becomes an offset vector in the direction of the side that corresponds to the x axis in the net, the y is in the direction of the height vector to the remaining vertex.
    Normalize the resulting point to project it to the unit sphere.
    Convert that to latitude and longitude.
    Get that pixel value from the source data for that lat/long.
    
It would be better to use the containing icosohedron rather than the contained icosohedron (so he centers of the faces touch the unit spehere instead of the vertices).
"""

import math
from PIL import Image, ImageDraw

SQRT3 = math.sqrt(3)

ctheta = math.cos(math.atan(0.5))
stheta = math.sin(math.atan(0.5))
IcosVertices = [(0, 0, 1),
                (ctheta * math.cos(math.radians(72*0)),    ctheta * math.sin(math.radians(72*0)),    stheta),
                (ctheta * math.cos(math.radians(72*1)),    ctheta * math.sin(math.radians(72*1)),    stheta),
                (ctheta * math.cos(math.radians(72*2)),    ctheta * math.sin(math.radians(72*2)),    stheta),
                (ctheta * math.cos(math.radians(72*3)),    ctheta * math.sin(math.radians(72*3)),    stheta),
                (ctheta * math.cos(math.radians(72*4)),    ctheta * math.sin(math.radians(72*4)),    stheta),
                (ctheta * math.cos(math.radians(36+72*0)), ctheta * math.sin(math.radians(36+72*0)), -stheta),
                (ctheta * math.cos(math.radians(36+72*1)), ctheta * math.sin(math.radians(36+72*1)), -stheta),
                (ctheta * math.cos(math.radians(36+72*2)), ctheta * math.sin(math.radians(36+72*2)), -stheta),
                (ctheta * math.cos(math.radians(36+72*3)), ctheta * math.sin(math.radians(36+72*3)), -stheta),
                (ctheta * math.cos(math.radians(36+72*4)), ctheta * math.sin(math.radians(36+72*4)), -stheta),
                (0, 0, -1),
                ]

IcosFaces = [(0,2,1), (0,3,2), (0,4,3), (0,5,4), (0,1,5),
             (1,2,6), (2,3,7), (3,4,8), (4,5,9), (5,1,10),
             (2,7,6), (3,8,7), (4,9,8), (5,10,9), (1,6,10),
             (6,7,11), (7,8,11), (8,9,11), (9,10,11), (10,6,11),
             ]


class IcosoProjection:
    """Do our terrible brute force icosohedral map."""
    def __init__(self, width=1000):
        """Arguments:
               width(int): width of the image
        """
        self.width = width
        self.height = int(math.ceil((SQRT3 / 2.0 * (width / 5.5)) * 3))
        self.initialize()
    
    def initialize(self):
        """Initialize the width dependent data."""
        halfSide = (self.width -1) / 11.0
        height = SQRT3 / 2.0 * (halfSide * 2.0)
        triangles = [
                     [(1, 0), (2, 1), (0, 1)],
                     [(3, 0), (4, 1), (2, 1)],
                     [(5, 0), (6, 1), (4, 1)],
                     [(7, 0), (8, 1), (6, 1)],
                     [(9, 0), (10, 1), (8, 1)],
                     [(0, 1), (2, 1), (1, 2)],
                     [(2, 1), (4, 1), (3, 2)],
                     [(4, 1), (6, 1), (5, 2)],
                     [(6, 1), (8, 1), (7, 2)],
                     [(8, 1), (10, 1), (9, 2)],
                     [(2, 1), (3, 2), (1, 2)],
                     [(4, 1), (5, 2), (3, 2)],
                     [(6, 1), (7, 2), (5, 2)],
                     [(8, 1), (9, 2), (7, 2)],
                     [(10, 1), (11, 2), (9, 2)],
                     [(1, 2), (3, 2), (2, 3)],
                     [(3, 2), (5, 2), (4, 3)],
                     [(5, 2), (7, 2), (6, 3)],
                     [(7, 2), (9, 2), (8, 3)],
                     [(9, 2), (11, 2), (10, 3)],
                     ]
        self.triangles = []
        for triangle in triangles:
            tri= []
            for x, y in triangle:
                tri.append((int(round(x * halfSide)), int(round(y * height))))
            self.triangles.append(tri)
            
    def overlayTriangles(self, image):
        """Draw the triangles.
        Arguments:
            image(PIL Image): the image to draw on, must match the width/height.
        """
        draw = ImageDraw.Draw(image)
        # the \ lines
        draw.line((self.triangles[0][2], self.triangles[15][2]), fill='black', width=1)
        draw.line((self.triangles[0][0], self.triangles[16][2]), fill='black', width=1)
        draw.line((self.triangles[1][0], self.triangles[17][2]), fill='black', width=1)
        draw.line((self.triangles[2][0], self.triangles[18][2]), fill='black', width=1)
        draw.line((self.triangles[3][0], self.triangles[19][2]), fill='black', width=1)
        draw.line((self.triangles[4][0], self.triangles[19][1]), fill='black', width=1)
        # the / lines
        draw.line((self.triangles[0][0], self.triangles[0][2]), fill='black', width=1)
        draw.line((self.triangles[1][0], self.triangles[15][0]), fill='black', width=1)
        draw.line((self.triangles[2][0], self.triangles[15][2]), fill='black', width=1)
        draw.line((self.triangles[3][0], self.triangles[16][2]), fill='black', width=1)
        draw.line((self.triangles[4][0], self.triangles[17][2]), fill='black', width=1)
        draw.line((self.triangles[4][1], self.triangles[18][2]), fill='black', width=1)
        draw.line((self.triangles[19][1], self.triangles[19][2]), fill='black', width=1)
        # the - lines
        draw.line((self.triangles[0][2], self.triangles[4][1]), fill='black', width=1)
        draw.line((self.triangles[15][0], self.triangles[19][1]), fill='black', width=1)
        del draw

    def createMap(self, dataSource):
        """Create the map.
        
        Arguments:
            dataSource(callable): when passed long, lat should return the color to plot as (R, G, B, A).
        Returns:
            PIL Image of the map.
        """
        image = Image.new("RGBA", (self.width, self.height))
        pix = image.load()
        for y in range(self.height):
            for x in range(self.width):
                ll = self.calcLatLong(x, y)
                if ll:
                    pix[x, y] = dataSource(ll[0], ll[1])
        return image
    
    def calcLatLong(self, x, y):
        """Get Lat, Long for a point on the map.
        
        Arguments:
            x(int): the x-coordinate of the pixel
            y(int): the y-coordinate of the pixel
        Returns:
            (long, lat) if (x, y) is inside the map portion of the image, None if not.
        """
        triangle = self.getContainingTriangle(x, y)
        if triangle is None:
            return None
        offset = self.triangleOffset(triangle, x, y)
        icosCoord = self.mapToIcosohedronTriangle(triangle, offset)
        return self.coordToLatLong(icosCoord)

    
    def getContainingTriangle(self, x, y):
        """Get the triangle the point is in."""
        def sameSide(p1, p2, a, b):
            """Are p1 and p2 on the same side of the line a, b?"""
            return 0 <= ((p1[1] - a[1]) * (b[0] - a[0]) - (p1[0] - a[0]) * (b[1] - a[1])) * ((p2[1] - a[1]) * (b[0] - a[0]) - (p2[0] - a[0]) * (b[1] - a[1]))
        
        for i in range(len(self.triangles)):
            a, b, c = self.triangles[i]
            if (x,y) == a or (x,y) == b or (x,y) == c or sameSide((x, y), a, b, c) and sameSide((x, y), b, a, c) and sameSide((x, y), c, a, b):
                return i
        return None
    
    def triangleOffset(self, triangle, x, y):
        """get offset from left corner."""
        if triangle / 5 in (0, 2):
            xoffset = (x - self.triangles[triangle][2][0]) / float(self.triangles[triangle][1][0] - self.triangles[triangle][2][0])
            yoffset = (y - self.triangles[triangle][2][1]) / float(self.triangles[triangle][0][1] - self.triangles[triangle][2][1])
        else:
            xoffset = (x - self.triangles[triangle][0][0]) / float(self.triangles[triangle][1][0] - self.triangles[triangle][0][0])
            yoffset = (y - self.triangles[triangle][0][1]) / float(self.triangles[triangle][2][1] - self.triangles[triangle][0][1])
        return xoffset, yoffset
    
    def mapToIcosohedronTriangle(self, triangle, offset):
        """Given the triangle and x,y offset on the 2D nest return the (x, y, z) for the corresponding point on the icosohedron face."""
        if triangle / 5 in (0, 2):
            offsetVertex = IcosVertices[IcosFaces[triangle][2]]
            offsetYVectorEnd = IcosVertices[IcosFaces[triangle][0]]

        else:
            offsetVertex = IcosVertices[IcosFaces[triangle][0]]
            offsetYVectorEnd = IcosVertices[IcosFaces[triangle][2]]
            
        offsetXVector = (IcosVertices[IcosFaces[triangle][1]][0] - offsetVertex[0],
                         IcosVertices[IcosFaces[triangle][1]][1] - offsetVertex[1],
                         IcosVertices[IcosFaces[triangle][1]][2] - offsetVertex[2])
        offsetYVector = (offsetYVectorEnd[0] - offsetVertex[0] - offsetXVector[0]/2.0,
                         offsetYVectorEnd[1] - offsetVertex[1] - offsetXVector[1]/2.0,
                         offsetYVectorEnd[2] - offsetVertex[2] - offsetXVector[2]/2.0)
        
        return (offsetVertex[0] + offset[0] * offsetXVector[0] + offset[1] * offsetYVector[0],
                offsetVertex[1] + offset[0] * offsetXVector[1] + offset[1] * offsetYVector[1],
                offsetVertex[2] + offset[0] * offsetXVector[2] + offset[1] * offsetYVector[2])


    def coordToLatLong(self, point):
        """Given a point (x, y, z) on the surface of a unit sphere centered at (0, 0, 0) return the latitude, longitude."""
        length = math.sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2])
        point = (point[0]/length, point[1]/length, point[2]/length)
        return math.asin(point[2]), math.atan2(point[1],  point[0]), 


class EquirectangularSource:
    """Using an equirectangular projection image provide a (lat, lon) interface to the points."""
    def __init__(self, image):
        self.image = image
        width, height = image.size
        self.widthStep = width / (2 * math.pi)
        self.heightStep = height / math.pi
        self.pixels = image.load()
        self.width = width
        self.height = height
        
    def latLongPixelValue(self, lat, lon):
        y = int(self.height / 2 - self.heightStep * lat)
        if y < 0:
            y = 0
        if y >= self.height:
            y = self.height - 1
        pi2 = math.pi + math.pi
        while lon < 0:
            lon = lon + pi2
        while lon > pi2:
            lon = lon- pi2
        x = int(self.widthStep * lon)
        if x < 0:
            x = 0
        if x >= self.width:
            x = self.width - 1
        return self.pixels[x, y]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="filename of source PNG image - equirectangular projection")
    parser.add_argument("output", help="filename to use for output PNG image")
    parser.add_argument("-w", "--width",  help="width of output image", type=int, default=1000)
    parser.add_argument("-g", "--grid",  help="draw triangle borders", action="store_true")

    args = parser.parse_args()
    source = EquirectangularSource(Image.open(args.input))         
    test = IcosoProjection(width=args.width)
    image = test.createMap(source.latLongPixelValue)
    if args.grid:
        test.overlayTriangles(image)
    image.save(args.output, "PNG")
