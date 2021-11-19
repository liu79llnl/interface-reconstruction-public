import vtk
import os

"""
#Write facet objects as vtk files
def writeFacets(facets, savename):
    try:
        os.mkdir("advection_vtk/{}".format(savename))
    except:
        pass

    g = open("advection_vtk/{}_all.visit".format(savename), "w")

    facetnum = 0
    for facet in facets:
        if facet.name == 'linear':
            line = vtk.vtkLineSource()
            line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
        elif facet.name == 'arc':
            arc = vtk.vtkArcSource()
            arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            arc.SetCenter(facet.center[0], facet.center[1], 0)
            arc.SetResolution(8)
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
            writer.SetInputConnection(arc.GetOutputPort())
            writer.Write()
            facetnum += 1
        elif facet.name == 'corner':
            if facet.radiusLeft is not None:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                arc.SetPoint2(facet.corner[0], facet.corner[1], 0)
                arc.SetCenter(facet.centerLeft[0], facet.centerLeft[1], 0)
                arc.SetResolution(8)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(arc.GetOutputPort())
            else:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
                line.SetPoint2(facet.corner[0], facet.corner[1], 0)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1
            if facet.radiusRight is not None:
                arc = vtk.vtkArcSource()
                arc.SetPoint1(facet.corner[0], facet.corner[1], 0)
                arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                arc.SetCenter(facet.centerRight[0], facet.centerRight[1], 0)
                arc.SetResolution(8)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(arc.GetOutputPort())
            else:
                line = vtk.vtkLineSource()
                line.SetPoint1(facet.corner[0], facet.corner[1], 0)
                line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
                writer = vtk.vtkXMLPolyDataWriter()
                writer.SetFileName("advection_vtk/{}/facet_{}.vtp".format(savename, facetnum))
                writer.SetInputConnection(line.GetOutputPort())
            writer.Write()
            facetnum += 1

    g.write("!NBLOCKS {}\n".format(facetnum))
    for i in range(facetnum):
        g.write("{}/facet_{}.vtp\n".format(savename, i))
    g.close()
"""

def writeFacets(facets, savename):

    vtkappend = vtk.vtkAppendPolyData()

    for facet in facets:
        if facet.name == 'linear':
            line = vtk.vtkLineSource()
            line.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            line.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            line.Update()
            vtkappend.AddInputData(line.GetOutput())
        elif facet.name == 'arc':
            arc = vtk.vtkArcSource()
            arc.SetPoint1(facet.pLeft[0], facet.pLeft[1], 0)
            arc.SetPoint2(facet.pRight[0], facet.pRight[1], 0)
            arc.SetCenter(facet.center[0], facet.center[1], 0)
            arc.SetResolution(8)
            arc.Update()
            vtkappend.AddInputData(arc.GetOutput())
    
    vtkappend.Update()
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName("advection_vtk/{}.vtp".format(savename))
    writer.SetInputConnection(vtkappend.GetOutputPort())
    writer.Update()
    writer.Write()