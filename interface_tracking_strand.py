import math

from mesh import QuadMesh, makeFineCartesianGrid, makeQuadGrid
from initialize_areas import initializeCircle, initializePoly, zalesak, xpluso
from interface_reconstruction import merge, makeFacets

from facet import getNormal, LinearFacet, ArcFacet
from strand import StrandPoly, Strand
from write_facets import writeFacets

import matplotlib.pyplot as plt
from geoms import getArea

#Hyperparameters
totalt = 8
dt = 0.0025

checkSize = 2

makeGapless = True

#Test settings
threshold = 1e-10

meshSize = 100
resolution = 32/100
meshType = 'cartesian' #quads #cartesian

#testType = 'vortex' #zalesak, xpluso, deformation
testType = 'vortex'

#----------------
#Code

#Initialize test settings
timesteps = int(totalt/dt)

if meshType == 'quads':
    opoints = makeQuadGrid(meshSize, resolution)
elif meshType == 'cartesian':
    opoints = makeFineCartesianGrid(meshSize, resolution)
mesh = QuadMesh(opoints)
opolys = mesh.polys

if testType == 'vortex':
    areas = initializeCircle(opolys, [50.01, 75.01], 15, threshold)
    velocity = lambda t, p : [-200*math.cos(math.pi*t/totalt)*(math.sin(math.pi*p[0]/100))**2 * math.sin(math.pi*p[1]/100) * math.cos(math.pi*p[1]/100),
                              200*math.cos(math.pi*t/totalt)*math.sin(math.pi*p[0]/100)*math.cos(math.pi*p[0]/100) * (math.sin(math.pi*p[1]/100))**2]
if testType == 'deformation':
    areas = initializeCircle(opolys, [50.01, 75.01], 15, threshold)
    velocity = lambda t, p: [100*math.cos(math.pi*t/totalt)*math.sin(4*math.pi*(p[0]+50)/100)*math.sin(4*math.pi*(p[1]+50)/100),
                             100*math.cos(math.pi*t/totalt)*math.cos(4*math.pi*(p[0]+50)/100)*math.cos(4*math.pi*(p[1]+50)/100)]
elif testType == 'zalesak':
    areas = zalesak(opolys, threshold)
    velocity = lambda t, p: [50-p[1], p[0]-50]
elif testType == 'xpluso':
    areas = xpluso(opolys, threshold)
    velocity = lambda t, p: [0.1, 0.1]

area_fractions = [[0 for _ in range(len(areas[0]))] for _ in range(len(areas))]
for x in range(len(areas)):
    for y in range(len(areas)):
        area_fractions[x][y] = areas[x][y] / getArea(opolys[x][y])
        if area_fractions[x][y] < threshold:
            area_fractions[x][y] = 0
        elif area_fractions[x][y] > 1-threshold:
            area_fractions[x][y] = 1

mesh.setAreas(area_fractions)
mesh.plotPartialAreas('advection_plt/areas_initial.png')

#------------------

#Interface reconstruction: compute initial interface
mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions = merge(opolys, areas)
predfacets, facetsunique = makeFacets(mergedpolyindices, mergedpolyinfos, mergedcoords, mergedareafractions, opolys, areas, threshold, makeGapless)


#for facetunique in facetsunique:
#    if facetunique is not None:
#        print(facetunique)

#Plot mesh
mesh.plotMesh('advection_vtk/quads_{}x{}.vtk'.format(int(meshSize*resolution), int(meshSize*resolution)))


#Plot initial setup
#writeFacets(facetsunique, 'simple_advect/facets_initial')

#------------------
#Advection only
"""
newfacetsunique = []
for facet in facetsunique:
    newfacet = ArcFacet(facet[1], facet[2], facet[3][0], facet[3][1])
    newfacetsunique.append(newfacet)
facetsunique = newfacetsunique
"""
"""
#Direct advection
t = 0
for timestep in range(timesteps):
    print("Direct advection: timestep {}".format(timestep))
    newfacetsunique = []
    for facet in facetsunique:
        newfacet = facet.advected(velocity, t, dt)
        newfacetsunique.append(newfacet)
    facetsunique = newfacetsunique
    writeFacets(facetsunique, 'simple_advect/facets_timestep_{}'.format(timestep))
    t += dt
"""

#------------------

#Convert to interface tracking format: StrandPolys
strandpolys = [[None] * len(predfacets[0]) for _ in range(len(predfacets))]
for x in range(len(strandpolys)):
    for y in range(len(strandpolys[0])):
        strandpolys[x][y] = StrandPoly(opolys[x][y])
        
        if x == 62 and y == 77:
            print(opolys[x][y])
        #Add facet info
        facet = predfacets[x][y]
        if facet is not None and facet[0] == 'arc':
            newfacet = ArcFacet(facet[1], facet[2], facet[3][0], facet[3][1])
            newstrand = Strand([newfacet])
            strandpolys[x][y].gatherStrands(newstrand)
            #print(strandpolys[x][y].strands)
            #print(strandpolys[x][y].strands_poly_lerps)
            strandpolys[x][y].updateStrands()
            strandpolys[x][y].zeroStrands()

            #print(strandpolys[x][y])


"""
#Attempt to initialize x+o facets
strandpolys = [[None] * len(opolys[0]) for _ in range(len(opolys))]
all_initial_facets = []
init_cross = [[7, 3], [13, 9], [19, 3], [23, 7], [17, 13], [23, 19], [19, 23], [13, 17], [7, 23], [3, 19], [9, 13], [3, 7]]
for i in range(len(init_cross)):
    all_initial_facets.append(LinearFacet(init_cross[i], init_cross[(i+1) % len(init_cross)]))
init_cross = [[9, 28], [17, 28], [17, 34], [23, 34], [23, 42], [17, 42], [17, 48], [9, 48], [9, 42], [3, 42], [3, 34], [9, 34]]
for i in range(len(init_cross)):
    all_initial_facets.append(LinearFacet(init_cross[i], init_cross[(i+1) % len(init_cross)]))
init_circle_points = [[48, 13], [38, 23], [28, 13], [38, 3]]
for i in range(4):
    all_initial_facets.append(ArcFacet([38, 13], 10, init_circle_points[i], init_circle_points[(i+1) % 4]))
init_circle_points = [[44, 13], [38, 7], [32, 13], [38, 19]]
for i in range(4):
    all_initial_facets.append(ArcFacet([38, 13], -6, init_circle_points[i], init_circle_points[(i+1) % 4]))
init_circle_points = [[41, 13], [38, 16], [35, 13], [38, 10]]
for i in range(4):
    all_initial_facets.append(ArcFacet([38, 13], 3, init_circle_points[i], init_circle_points[(i+1) % 4]))
init_circle_points = [[48, 38], [38, 48], [28, 38], [38, 28]]
for i in range(4):
    all_initial_facets.append(ArcFacet([38, 38], 10, init_circle_points[i], init_circle_points[(i+1) % 4]))
init_circle_points = [[44, 38], [38, 32], [32, 38], [38, 44]]
for i in range(4):
    all_initial_facets.append(ArcFacet([38, 38], -6, init_circle_points[i], init_circle_points[(i+1) % 4]))

for x in range(len(strandpolys)):
    for y in range(len(strandpolys[0])):
        strandpolys[x][y] = StrandPoly(opolys[x][y])
        for all_initial_facet in all_initial_facets:
            newstrand = Strand([all_initial_facet])
            strandpolys[x][y].gatherStrands(newstrand)
            strandpolys[x][y].updateStrands()
            strandpolys[x][y].zeroStrands()
"""

initial_facets = []
for x in range(len(strandpolys)):
    for y in range(len(strandpolys[0])):
        for strand in strandpolys[x][y].strands_stored:
            for facet in strand.facets:
                initial_facets.append(facet)
writeFacets(initial_facets, 'facets_initial')

"""
facetclass = [[None] * len(predfacets[0]) for _ in range(len(predfacets))]
for x in range(len(predfacets)):
    for y in range(len(predfacets[0])):
        facet = predfacets[x][y]
        if facet[0] == 'linear':
            newfacet = LinearFacet(facet[1], facet[2])
        elif facet[0] == 'arc':
            newfacet = ArcFacet(facet[1], facet[2], facet[3][0], facet[3][1])
        elif facet[0] == 'corner':
            newfacet = CornerFacet(None, None, None, None, facet[1][0], facet[1][1], facet[1][2])
        elif facet[0] == 'curvedcorner':
            newfacet = CornerFacet(facet[1], facet[2], facet[3], facet[4], facet[5][0], facet[5][1], facet[5][2])
        facetclass[x][y] = newfacet
"""

print('Advection')
#timesteps = 520
t = 0
for timestep in range(timesteps):
    print("Advection: timestep {}".format(timestep))
    for x in range(len(strandpolys)):
        for y in range(len(strandpolys[0])):
            #Compute advected strands
            strands_stored = strandpolys[x][y].strands_stored
            advected_strands = []
            for strand in strands_stored:
                advected_strands.append(strand.advectedStrand(velocity, t, dt))
            #Neighbors
            for checkSize_x in range(-checkSize, checkSize+1):
                for checkSize_y in range(-checkSize, checkSize+1):
                    check_x = x+checkSize_x
                    check_y = y+checkSize_y
                    if check_x == 62 and check_y == 77 and timestep == 521:
                        for advected_strand in advected_strands:
                            print(advected_strand)
                    if check_x >= 0 and check_x < len(strandpolys) and check_y >= 0 and check_y < len(strandpolys[0]):
                        #Valid neighbor of [x, y]
                        for advected_strand in advected_strands:
                            strandpolys[check_x][check_y].gatherStrands(advected_strand)
                            
    """
    #Print premerged strands
    print("Premerged timestep {}".format(timestep))
    for x in range(len(strandpolys)):
        for y in range(len(strandpolys[0])):
            for strand in strandpolys[x][y].strands:
                print(strandpolys[x][y].points)
                print(strand)
    """
    """
    #Plot premerged strands
    if timestep % 25 == 24 or (timestep < 5 and timestep > -1):
        premergedFacets = []
        for x in range(len(strandpolys)):
            for y in range(len(strandpolys[0])):
                for strand in strandpolys[x][y].strands:
                    for facet in strand.facets:
                        premergedFacets.append(facet)
        writeFacets(premergedFacets, 'simple_advect/facets_premerged_{}'.format(timestep))
    """
    print("Begin merging")
    #Merge strands
    for x in range(len(strandpolys)):
        for y in range(len(strandpolys[0])):
            if x == 62 and y == 77 and timestep == 521:
                for advected_strand in strandpolys[x][y].merged_strands:
                    print("Merged")
                    print(advected_strand)
                for advected_strand in strandpolys[x][y].unmerged_strands:
                    print("unmerged")
                    print(advected_strand)
            if timestep == timesteps - 1:
                strandpolys[x][y].updateStrands(merge_all=True)
            else:
                strandpolys[x][y].updateStrands(merge_all=False)
            #Zero
            strandpolys[x][y].zeroStrands()
            if x == 62 and y == 77 and timestep == 521:
                print("Merged")
                for advected_strand in strandpolys[x][y].strands_stored:
                    print(advected_strand)


    #Plot facets
    if timestep % 1 == 0 or (timestep < 5 and timestep > -1):
        mergedFacets = []
        for x in range(len(strandpolys)):
            for y in range(len(strandpolys[0])):
                num_postmerged_facets = 0
                for strand in strandpolys[x][y].strands_stored:
                    for facet in strand.facets:
                        mergedFacets.append(facet)
                        num_postmerged_facets += 1
                if num_postmerged_facets > 1:
                    print(num_postmerged_facets)
        print(len(mergedFacets))

        writeFacets(mergedFacets, 'facets_merged_{}'.format(timestep))

    """
    #Print strands
    print("Postmerged timestep {}".format(timestep))
    if timestep % 1 == 0:
        for x in range(len(strandpolys)):
            for y in range(len(strandpolys[0])):
                if len(strandpolys[x][y].strands_stored) > 0:
                    print(strandpolys[x][y].points)
                    for strand in strandpolys[x][y].strands_stored:
                        print(strand)
                    #print(strandpolys[x][y].strands_poly_lerps)
                    #strandpolys[x][y].updateStrands()
                    #print(strandpolys[x][y])
                continue
    """

    t += dt

curvature_errors = []
for x in range(len(strandpolys)):
    for y in range(len(strandpolys[0])):
        if len(strandpolys[x][y].strands_stored) > 0:
            for strand in strandpolys[x][y].strands_stored:
                for facet in strand.facets:
                    curvature_errors.append(math.log10(abs(facet.curvature - 1/15)))

plt.hist(curvature_errors)
plt.xlabel('log curvature error')
plt.ylabel('Frequency')
plt.title('Vortex test, {}x{} mesh, curvature errors'.format(int(meshSize*resolution), int(meshSize*resolution)))
plt.savefig('advection_plt/curvature_errors_{}x{}.png'.format(int(meshSize*resolution), int(meshSize*resolution)))
plt.clf()

#Compute areas
area_errors = []
interface_tracking_areas = [[None] * len(predfacets[0]) for _ in range(len(predfacets))]
for x in range(len(interface_tracking_areas)):
    for y in range(len(interface_tracking_areas[0])):
        interface_tracking_areas[x][y] = strandpolys[x][y].getArea()
        #compare ending area to beginning area
        if areas[x][y] != 1 and areas[x][y] != 0:
            print("Original area: {}".format(areas[x][y]))
            print("New area: {}".format(interface_tracking_areas[x][y]))
            area_errors.append(math.log10(abs(areas[x][y] - interface_tracking_areas[x][y])))

plt.hist(area_errors)
plt.xlabel('log area error')
plt.ylabel('Frequency')
plt.title('Vortex test, {}x{} mesh, area errors'.format(int(meshSize*resolution), int(meshSize*resolution)))
plt.savefig('advection_plt/area_errors_{}x{}.png'.format(int(meshSize*resolution), int(meshSize*resolution)))
plt.clf()



#------------------
#Interface tracking

#convert from format into strands per cell

#new formatting for strands in cell: sorted per initial endpoint, etc

#for each timestep:
#   advect the strands
#   separate the strands
#   remerge the strands


#advect the strands
#for each strand, for each facet in strand
#advect facet

#separate the strands
#???

#remerge the strands
#get areas of all arcs
#run arc thing