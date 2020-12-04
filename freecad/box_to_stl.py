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

Mesh.export([pf1, pf2], "bar.stl")
