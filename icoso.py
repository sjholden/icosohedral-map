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
If the point is inside one of the triangles (we use the simple check if the point is on the same side of the line defined by each vertex pair as the remaining vertex):
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

ctheta = math.cos(math.atan(2.0))
stheta = math.sin(math.atan(0.5))
# The vertices of our isocohedron.
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
# The faces of our isocohedron, each entry is the vertices of the triangle as indexes into IcosVertices.
IcosFaces = [(0,2,1), (0,3,2), (0,4,3), (0,5,4), (0,1,5),
             (1,2,6), (2,3,7), (3,4,8), (4,5,9), (5,1,10),
             (2,7,6), (3,8,7), (4,9,8), (5,10,9), (1,6,10),
             (6,7,11), (7,8,11), (8,9,11), (9,10,11), (10,6,11),
             ]

# The triangle mesh - the triangles must map to the corresponding Icosfaces - index position is used to
# match them up. Laid out on a 12x4 grid with the vertices on integer coordinates - we stretch to be
# equaliateral triangles before doing any mapping.
TriangleMesh = [[(1, 0), (2, 1), (0, 1)],
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

# We walk the output image to populate, moving left to right and top to bottom, given our triangle mesh
# a step one pixel right (wrapping to the left and down one) can only reach a small number of
# triangles from the last triangle we visited, this list encodes that - so we check a few triangles
# instead of 20 at each step. There are more than just 2 entries for the middle rows since the rounding
# to pixels can jump across a corner in some cases.
# All indexes into TriangleMesh.
TrianglesNext = [(0,1),
                 (1,2),
                 (2,3),
                 (3,4),
                 (4,0,5),
                 (5,10,6),
                 (6,11,7),
                 (7,12,8),
                 (8,13,9),
                 (9,14,5,15),
                 (10,6,11),
                 (11,7,12),
                 (12,8,13),
                 (13,9,14),
                 (14,5,15),
                 (15,16),
                 (16,17),
                 (17,18),
                 (18,19),
                 (19,15),
                ]

# Our triangle mesh is this, labels are the indexes into TriangleMesh
# Note that (index // 5) % 2 distinguishes up from down triangles.
#   /\    /\    /\    /\    /\
#  / 0\  / 1\  / 2\  / 3\  / 4\
# /____\/____\/____\/____\/____\
# \  5 /\  6 /\  7 /\  8 /\  9 /\
#  \  /10\  /11\  /12\  /13\  /14\
#   \/____\/____\/____\/____\/____\
#    \ 15 /\ 16 /\ 17 /\ 18 /\ 19 /
#     \  /  \  /  \  /  \  /  \  /
#      \/    \/    \/    \/    \/




class IcosoProjection:
    """Do our terrible brute force icosohedral map."""
    def __init__(self, width=1000):
        """Arguments:
               width(int): width of the image
        """
        self.width = width
        self.height = int(math.ceil((SQRT3 / 2.0 * (width / 5.5)) * 3))
        self.initialize()
        self.currentTriangle = 0

    def initialize(self):
        """Initialize the width dependent data."""
        halfSide = (self.width -1) / 11.0
        height = SQRT3 / 2.0 * (halfSide * 2.0)
        self.triangles = []
        for triangle in TriangleMesh:
            tri = []
            for x, y in triangle:
                tri.append((x * halfSide, y * height))
            self.triangles.append(tri)

    def overlayTriangles(self, image):
        """Draw the triangles.
        Arguments:
            image(PIL Image): the image to draw on, must match the width/height.
        """
        draw = ImageDraw.Draw(image)
        # the \ lines
        draw.line((self.triangles[0][2], self.triangles[15][2]), fill='black', width=2)
        draw.line((self.triangles[0][0], self.triangles[16][2]), fill='black', width=2)
        draw.line((self.triangles[1][0], self.triangles[17][2]), fill='black', width=2)
        draw.line((self.triangles[2][0], self.triangles[18][2]), fill='black', width=2)
        draw.line((self.triangles[3][0], self.triangles[19][2]), fill='black', width=2)
        draw.line((self.triangles[4][0], self.triangles[19][1]), fill='black', width=2)
        # the / lines
        draw.line((self.triangles[0][0], self.triangles[0][2]), fill='black', width=2)
        draw.line((self.triangles[1][0], self.triangles[15][0]), fill='black', width=2)
        draw.line((self.triangles[2][0], self.triangles[15][2]), fill='black', width=2)
        draw.line((self.triangles[3][0], self.triangles[16][2]), fill='black', width=2)
        draw.line((self.triangles[4][0], self.triangles[17][2]), fill='black', width=2)
        draw.line((self.triangles[4][1], self.triangles[18][2]), fill='black', width=2)
        draw.line((self.triangles[19][1], self.triangles[19][2]), fill='black', width=2)
        # the - lines
        draw.line((self.triangles[0][2], self.triangles[4][1]), fill='black', width=2)
        draw.line((self.triangles[15][0], self.triangles[19][1]), fill='black', width=2)
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
        self.currentTriangle = 0
        # local variables instead of self lookups in the loop just for speed
        getContainingTriangle = self.getContainingTriangle
        triangles = self.triangles
        for y in range(self.height):
            for x in range(self.width):
                # Find the triangle we are in
                triangle = getContainingTriangle(x, y)
                if triangle is None:
                    continue
                # Find the offset within that triangle from the left corner - up and down triangles are different
                # Find the corresponding offset in the icosohedron face, in 3 dimensions
                if triangle // 5 in (0, 2):
                    xoffset = (x - triangles[triangle][2][0]) / (triangles[triangle][1][0] - triangles[triangle][2][0])
                    yoffset = (y - triangles[triangle][2][1]) / (triangles[triangle][0][1] - triangles[triangle][2][1])
                    offsetVertex = IcosVertices[IcosFaces[triangle][2]]
                    offsetYVectorEnd = IcosVertices[IcosFaces[triangle][0]]
                else:
                    xoffset = (x - triangles[triangle][0][0]) / (triangles[triangle][1][0] - triangles[triangle][0][0])
                    yoffset = (y - triangles[triangle][0][1]) / (triangles[triangle][2][1] - triangles[triangle][0][1])
                    offsetVertex = IcosVertices[IcosFaces[triangle][0]]
                    offsetYVectorEnd = IcosVertices[IcosFaces[triangle][2]]
                offsetXVector = (IcosVertices[IcosFaces[triangle][1]][0] - offsetVertex[0],
                                 IcosVertices[IcosFaces[triangle][1]][1] - offsetVertex[1],
                                 IcosVertices[IcosFaces[triangle][1]][2] - offsetVertex[2])
                offsetYVector = (offsetYVectorEnd[0] - offsetVertex[0] - offsetXVector[0]/2.0,
                                 offsetYVectorEnd[1] - offsetVertex[1] - offsetXVector[1]/2.0,
                                 offsetYVectorEnd[2] - offsetVertex[2] - offsetXVector[2]/2.0)
                # Convert those offsets to five the coordinates on the icosohedron
                iCoord = (offsetVertex[0] + xoffset * offsetXVector[0] + yoffset * offsetYVector[0],
                          offsetVertex[1] + xoffset * offsetXVector[1] + yoffset * offsetYVector[1],
                          offsetVertex[2] + xoffset * offsetXVector[2] + yoffset * offsetYVector[2])
                # Convert that point to a latitide, longitude 
                length = math.sqrt(iCoord[0]*iCoord[0] + iCoord[1]*iCoord[1] + iCoord[2]*iCoord[2])
                point = (iCoord[0]/length, iCoord[1]/length, iCoord[2]/length)
                lat = math.asin(point[2])
                lon = math.atan2(point[1], point[0])
                pix[x, y] = dataSource(lat, lon)
        return image

    def getContainingTriangle(self, x, y):
        """Get the triangle the point is in."""
        def sameSide(p1, p2, a, b):
            """Are p1 and p2 on the same side of the line a, b?"""
            return 0 <= ((p1[1] - a[1]) * (b[0] - a[0]) - (p1[0] - a[0]) * (b[1] - a[1])) * ((p2[1] - a[1]) * (b[0] - a[0]) - (p2[0] - a[0]) * (b[1] - a[1]))
        for i in TrianglesNext[self.currentTriangle]:
            a, b, c = self.triangles[i]
            if (x,y) == a or (x,y) == b or (x,y) == c or sameSide((x, y), a, b, c) and sameSide((x, y), b, a, c) and sameSide((x, y), c, a, b):
                self.currentTriangle = i
                return i
        return None

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
        lon = lon % (2 * math.pi)
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
    source = EquirectangularSource(Image.open(args.input).convert('RGBA'))
    test = IcosoProjection(width=args.width)
    image = test.createMap(source.latLongPixelValue)
    if args.grid:
        test.overlayTriangles(image)
    image.save(args.output, "PNG")
