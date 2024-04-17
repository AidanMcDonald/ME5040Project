"""
beamExample.py

Reproduce the cantilever beam example from the
Appendix of the Getting Started with
ABAQUS Manual.
"""

from abaqus import *
from abaqusConstants import *

# Create a model.

pvModel = mdb.Model(name='PressureVessel')

# Create a new viewport in which to display the model
# and the results of the analysis.

myViewport = session.Viewport(name='Pressure Vessel',
    origin=(20, 20), width=150, height=120)
    
#-----------------------------------------------------

import part

# Create a sketch for the base feature.

pvRotateSketch = pvModel.Sketch(name='pvRotateSketch',
    sheetSize=250.)

pvRotateSketch.Line(point1=(20, 40), point2=(20, -40))
pvRotateSketch.Line(point1=(20, 40), point2=(0, 40))
pvRotateSketch.ArcByCenterEnds(center=(0, -40), point1=(20, -40), point2=(0,-60))
pvRotateSketch.assignCenterline(pvRotateSketch.Line(point1=(0, 40), point2=(0, -40)))

# Create a three-dimensional, deformable part.

pvPart = pvModel.Part(name='PressureVessel', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)

# Create the part's base feature by extruding the sketch 
# through a distance of 25.0.

pvPart.BaseSolidRevolve(sketch=pvSketch, angle=360)


#-----------------------------------------------------

import material

# Create a material.

mySteel = myModel.Material(name='Steel')

# Create the elastic properties: youngsModulus is 209.E3
# and poissonsRatio is 0.3

elasticProperties = (209.E3, 0.3)
mySteel.Elastic(table=(elasticProperties, ) )

#-------------------------------------------------------

import section

# Create the solid section.

mySection = myModel.HomogeneousSolidSection(name='pvSection',
    material='Steel', thickness=1.0)

# Assign the section to the region. The region refers 
# to the single cell in this model.

region = (pvPart.cells,)
pvPart.SectionAssignment(region=region,
    sectionName='pvSection')

#-------------------------------------------------------


import assembly

# Create a part instance.

myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name='pvInstance',
    part=pvPart, dependent=OFF)

#-------------------------------------------------------

import step

# Create a step. The time period of the static step is 1.0, 
# and the initial incrementation is 0.1; the step is created
# after the initial step. 

myModel.StaticStep(name='pressureStep', previous='Initial',
    timePeriod=1.0, initialInc=0.1,
    description='Pressurize it.')

#-------------------------------------------------------

import load

# Create a pressure load on the top face of the beam.

topSurface = [(face, SIDE2) for face in myInstance.faces]
myModel.Pressure(name='Pressure', createStepName='pressureStep',
    region=topSurface, magnitude=0.5)

#-------------------------------------------------------

import mesh

# Assign an element type to the part instance.

region = (myInstance.cells,)
elemType = mesh.ElemType(elemCode=C3D8I, elemLibrary=STANDARD)
myAssembly.setElementType(regions=region, elemTypes=(elemType,))

# Seed the part instance.

myAssembly.seedPartInstance(regions=(myInstance,), size=10.0)

# Mesh the part instance.

myAssembly.generateMesh(regions=(myInstance,))

# Display the meshed beam.

myViewport.assemblyDisplay.setValues(mesh=ON)
myViewport.assemblyDisplay.meshOptions.setValues(meshTechnique=ON)
myViewport.setValues(displayedObject=myAssembly)

#-------------------------------------------------------
'''
# Set unique material properties for each cell
for e in myAssembly.instances['pvInstance'].elements:

    # Create material based on element's location
    steelTemp = myModel.Material(name='Steel'+str(e.label))

    # Create the elastic properties: youngsModulus is 209.E3
    # and poissonsRatio is 0.3

    elasticProperties = (209.E3, 0.3)
    steelTemp.Elastic(table=(elasticProperties, ) )


    # Create section with region equal to the element region
    
    mySection = myModel.HomogeneousSolidSection(name='pvSection'+str(e.label),
        material='Steel'+str(e.label), thickness=1.0)

    region = regionToolset.Region(elements=e)
    pvPart.SectionAssignment(region=region,
        sectionName='pvSection'+str(e.label))
'''

#-------------------------------------------------------

import job

# Create an analysis job for the model and submit it.

jobName = 'pressure_vessel'
myJob = mdb.Job(name=jobName, model='PressureVessel',
    description='pressure vessel project')

# Wait for the job to complete.

myJob.submit()
myJob.waitForCompletion()

#-------------------------------------------------------

import visualization

# Open the output database and display a
# default contour plot.

myOdb = visualization.openOdb(path=jobName + '.odb')
myViewport.setValues(displayedObject=myOdb)
myViewport.odbDisplay.setPlotMode(CONTOUR)
