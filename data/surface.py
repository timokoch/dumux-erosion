import gmsh
import sys

gmsh.initialize(sys.argv)

# this example shows how alternative non-plane geometries can be used in the
# built-in kernel; here using a parametric surface

SURFACES = {
    "surface.msh": "0.2*Sin(Pi*u)*Sin(Pi*v)",
    "surface_2.msh": "0.25*Sin(Pi*u)*Sin(Pi*v)*Cos(0.33*Pi*u*v)*Tanh(4*Cos(u^2+v^2))",
    "surface_3.msh": "0.2*Sin(sqrt((5*u)^2+(5*v)^2))",
    "surface_4.msh": "0.1*Sin(sqrt((5*(u-2.5))^2+(5*v)^2)) + 0.1*Sin(sqrt((5*(u+2.5))^2+(5*v)^2))"
}

surface = "surface.msh"

g=gmsh.model.geo.addGeometry("ParametricSurface", strings=["u", "v", SURFACES[surface]])

p1=gmsh.model.geo.addPointOnGeometry(g, -5,-5,0,tag=1)
p2=gmsh.model.geo.addPointOnGeometry(g, 5,-5,0,tag=2)
p3=gmsh.model.geo.addPointOnGeometry(g, 5,5,0,tag=3)
p4=gmsh.model.geo.addPointOnGeometry(g, -5,5,0,tag=4)

l1=gmsh.model.geo.addLine(p1, p2)
l2=gmsh.model.geo.addLine(p2, p3)
l3=gmsh.model.geo.addLine(p3, p4)
l4=gmsh.model.geo.addLine(p4, p1)

cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.05)
gmsh.model.mesh.generate(2)
gmsh.write(surface)
gmsh.finalize()
