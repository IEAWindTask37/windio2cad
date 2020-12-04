import FreeCAD
import Part
import Mesh

doc = App.newDocument("Doc")

shape1 = Part.makeSphere(1)
pf1 = doc.addObject("Part::Feature", "myShape1")
pf1.Shape = shape1

shape2 = Part.makeSphere(1)
shape2.translate((10, 10, 10))
pf2 = doc.addObject("Part::Feature", "myShape2")
pf2.Shape = shape2

m = FreeCAD.Matrix()
m.move((2, 2, 2))
shape3 = Part.makeSphere(1)
shape3.transformShape(m)
pf3 = doc.addObject("Part::Feature", "myShape3")
pf3.Shape = shape3

fuse2 = shape1.fuse(shape2)
pf_fuse2 = doc.addObject("Part::Feature", "myFuse4")
pf_fuse2.Shape = fuse2

Mesh.export([pf_fuse2], "foo.stl")