import numpy as np
import math
import vtk
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

class QuadMesh:

    def __init__(self, points, areas=None):
        #Points = [[x, y]]
        self.points = points
        #Compute max x and y
        self.max_x = -float('inf')
        self.max_y = -float('inf')
        for x in range(len(self.points)):
            for y in range(len(self.points[0])):
                if points[x][y][0] > self.max_x:
                    self.max_x = points[x][y][0]
                if points[x][y][1] > self.max_y:
                    self.max_y = points[x][y][1]

        self.polys = [[None] * (len(points[0])-1) for _ in range(len(points)-1)]
        self.patches = []
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                poly = [points[x][y], points[x+1][y], points[x+1][y+1], points[x][y+1]]
                self.polys[x][y] = poly
                patch = Polygon(np.array(poly), True)
                self.patches.append(patch)
        
        if areas is None:
            self.areas = [[None] * (len(self.polys[0])) for _ in range(len(self.polys))]
            self.patchareas = np.array([])
            self.patchpartialareas = np.array([])
            self.patchinitialareas = np.array([])
        else:
            self.initializeAreas(areas)

    def plotMesh(self, path):
        sgrid = vtk.vtkStructuredGrid()
        sgrid.SetDimensions([len(self.points), len(self.points[0]), 1])
        vtkpoints = vtk.vtkPoints()
        vtkpoints.Allocate(len(self.points)*len(self.points[0]))
        counter = 0
        for x in range(len(self.points)):
            for y in range(len(self.points[0])):
                vtkpoints.InsertPoint(counter, [self.points[x][y][0], self.points[x][y][1], 0])
                counter += 1
        sgrid.SetPoints(vtkpoints)
        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName(path)
        writer.SetInputData(sgrid)
        writer.Write()

    def initializeAreas(self, areas):
        self.initialareas = areas
        patchinitialareas = []

        #Flatten array
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                patchinitialareas.append(areas[x][y])

        self.patchinitialareas = np.array(patchinitialareas)
        self.setAreas(areas)

    def setAreas(self, areas):
        self.areas = areas
        patchareas = []
        patchpartialareas = []

        #Flatten array
        for x in range(len(self.polys)):
            for y in range(len(self.polys[0])):
                a = areas[x][y]
                patchareas.append(a)
                patchpartialareas.append(math.ceil(a - math.floor(a)))

        self.patchareas = np.array(patchareas)
        self.patchpartialareas = np.array(patchpartialareas)

    def plotAreas(self, path):
        patchcollection = PatchCollection(self.patches, cmap='jet') #jet
        patchcollection.set_array(self.patchareas)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(0, self.max_x)
        ax.set_ylim(0, self.max_y)
        ax.add_collection(patchcollection)
        plt.savefig(path, dpi=199)
        plt.close()
        plt.clf()

    def plotPartialAreas(self, path):
        patchcollection = PatchCollection(self.patches, cmap='jet') #jet
        patchcollection.set_array(self.patchpartialareas)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(0, self.max_x)
        ax.set_ylim(0, self.max_y)
        ax.add_collection(patchcollection)
        plt.savefig(path, dpi=199)
        plt.close()
        plt.clf()

    def plotInitialAreaCompare(self, path):
        patchcollection = PatchCollection(self.patches, cmap='jet') #jet
        patchcollection.set_array(np.abs(self.patchareas-self.patchinitialareas))
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(0, self.max_x)
        ax.set_ylim(0, self.max_y)
        ax.add_collection(patchcollection)
        plt.savefig(path, dpi=199)
        plt.close()
        plt.clf()

#Generate cartesian mesh, gridSize*resolution x gridSize*resolution grid
def makeFineCartesianGrid(gridSize, resolution):
    print("Making quad mesh")
    points = [[0] * (int(gridSize*resolution)+1) for _ in range(int(gridSize*resolution)+1)]
    for x in range(len(points)):
        for y in range(len(points)):
            points[x][y] = [x/resolution, y/resolution]
    print("Done")
    return points

#Generate quad mesh
def makeQuadGrid(gridSize, resolution, wiggle=0.25):
    rng = np.random.RandomState(0)
    print("Making quad mesh")
    points = [[0] * (int(gridSize*resolution)+1) for _ in range(int(gridSize*resolution)+1)]
    for x in range(len(points)):
        for y in range(len(points)):
            points[x][y] = [(x + wiggle*rng.rand())/resolution, (y + wiggle*rng.rand())/resolution]
    print("Done")
    return points

#Generate concave mesh
def makeConcaveGrid(gridSize, wiggle):
    print("Making quad mesh")
    points = [[0] * (gridSize+1) for _ in range(gridSize+1)]
    for x in range(gridSize+1):
        for y in range(gridSize+1):
            if (x+y) % 2 == 1:
                points[x][y] = [x-wiggle, y-wiggle]
            else:
                points[x][y] = [x, y]
    print("Done")
    return points