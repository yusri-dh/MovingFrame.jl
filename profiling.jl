using GeometryBasics;
using FileIO;
using LinearAlgebra;
using GaussianProcesses;
using MovingFrame;
import Makie;
using OwnTime;
using Profile;
const mv = MovingFrame
##11
meshfile = Base.download("https://github.com/libigl/libigl-tutorial-data/raw/master/bumpy-cube.obj")
shape = FileIO.load(FileIO.File{FileIO.format"OBJ"}(meshfile))
verts = decompose(Point3f0, shape)
faces = decompose(TriangleFace{Int}, shape)
vmat = MovingFrame.array_of_array_to_matrix(verts)
fmat = MovingFrame.array_of_array_to_matrix(faces)
##
mv.surface_area2(vmat,fmat)
mv.surface_area(vmat,fmat)
##
mv.mean_curvature_flow(verts,faces)
mv.mean_curvature_flow(vmat,fmat)

