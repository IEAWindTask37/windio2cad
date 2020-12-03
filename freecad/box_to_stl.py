import FreeCAD
import Part
import Mesh

shape = Part.makeSphere(1)
doc = App.newDocument("Doc")
pf = doc.addObject("Part::Feature", "myShape")
pf.Shape = shape
Mesh.export([pf], "bar.stl")
