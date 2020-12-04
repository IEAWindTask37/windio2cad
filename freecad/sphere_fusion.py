import FreeCAD
import Part
import Mesh

doc = App.newDocument("Doc")

radii = [1.1, 1.2, 1.3, 1.4, 1.5, 1.4, 1.3, 1.2, 1.1]

union = []
part_features = []
for i, r in enumerate(radii):
    shape_name = f"shape_{i}"
    shape = Part.makeSphere(r)
    shape.translate((0, 0, i))
    part_feature = doc.addObject("Part::Feature", shape_name)
    part_feature.Shape = shape
    union.append(shape)
    part_features.append(part_feature)

Mesh.export(part_features, "foo.stl")
