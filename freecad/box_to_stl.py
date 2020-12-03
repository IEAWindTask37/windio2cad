import Part

FreeCAD.ActiveDocument = App.newDocument("myDocument")
myPart = FreeCAD.ActiveDocument.addObject("Part::Feature", "myPart")
cube = Part.makeBox(2, 2, 2)
myPart.Shape = cube
myPart.Shape.exportStl("foo.stl")
