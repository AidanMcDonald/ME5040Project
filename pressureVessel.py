"""
pressureVessel.py

Model a simple reactor pressure vessel with continuously
changing material properties.
"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import inspect
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

def printMethods(object):
    # Abaqus has terrible documentation so this function is useful for debugging
    object_methods = [method_name for method_name in dir(object)]
    print("object_methods: " + str(object_methods))

def getElasticProperties(element):
    elasticProperties = (209.E3, 0.3)
    return elasticProperties

# Create a model.

pvModel = mdb.Model(name='PressureVessel')

# Create a new viewport in which to display the model
# and the results of the analysis.

myViewport = session.Viewport(name='Pressure Vessel',
    origin=(20, 20), width=150, height=120)
    
#-----------------------------------------------------

import part

# Create a sketch for the base feature.

pvRotateSketch = pvModel.ConstrainedSketch(name='pvRotateSketch',
    sheetSize=100.)


pvRotateSketch.Line(point1=(0, 40), point2=(20, 40))
pvRotateSketch.Line(point1=(20, 40), point2=(20, -40))
pvRotateSketch.ArcByCenterEnds(center=(0, -40), point1=(0,-60), point2=(20, -40))
pvRotateSketch.Line(point1=(0, -60), point2=(0, -55))
pvRotateSketch.ArcByCenterEnds(center=(0, -40), point1=(0,-55), point2=(15,-40))
pvRotateSketch.Line(point1=(15, -40), point2=(15, 35))
pvRotateSketch.Line(point1=(15, 35), point2=(0, 35))
pvRotateSketch.Line(point1=(0, 35), point2=(0, 40))
pvRotateSketch.assignCenterline(pvRotateSketch.ConstructionLine(point1=(0,0), point2=(0,1)))

sketch = mdb.models["PressureVessel"].convertAllSketches()

# Create a three-dimensional, deformable part.

pvPart = pvModel.Part(name='PressureVessel', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)

# Create the part's base feature by extruding the sketch 
# through a distance of 25.0.

pvPart.BaseSolidRevolve(sketch=pvRotateSketch, angle=360)


#-----------------------------------------------------

import material

# Create a material.

mySteel = pvModel.Material(name='Steel')

# Create the elastic properties: youngsModulus is 209.E3
# and poissonsRatio is 0.3

elasticProperties = (209.E3, 0.3)
mySteel.Elastic(table=(elasticProperties, ) )

#-------------------------------------------------------

import section

# Create the solid section.

mySection = pvModel.HomogeneousSolidSection(name='pvSection',
    material='Steel', thickness=1.0)

# Assign the section to the region. The region refers 
# to the single cell in this model.

region = (pvPart.cells,)
pvPart.SectionAssignment(region=region,
    sectionName='pvSection')

#-------------------------------------------------------


import assembly

# Create a part instance.

myAssembly = pvModel.rootAssembly
myInstance = myAssembly.Instance(name='pvInstance',
    part=pvPart, dependent=OFF)

#-------------------------------------------------------

import step

# Create a step. The time period of the static step is 1.0, 
# and the initial incrementation is 0.1; the step is created
# after the initial step. 

pvModel.StaticStep(name='pressureStep', previous='Initial',
    timePeriod=1.0, initialInc=0.1,
    description='Pressurize it.')

#-------------------------------------------------------

import load

# Create a pressure load on the top face of the beam.
insideFacePoint = (15,0,0)
insideFace = myInstance.faces.findAt((insideFacePoint,) )

insideTopFacePoint = (0,35,0)
insideTopFace = myInstance.faces.findAt((insideTopFacePoint,) )

insideBottomFacePoint = (0,-55,0)
insideBottomFace = myInstance.faces.findAt((insideBottomFacePoint,) )

insideSurfaces = ((insideFace, SIDE1), (insideTopFace, SIDE1), (insideBottomFace, SIDE1))

pvModel.Pressure(name='Pressure', createStepName='pressureStep',
    region=insideSurfaces, magnitude=0.5)

#-------------------------------------------------------

import mesh

# Assign an element type to the part instance.

region = [cell for cell in myInstance.cells]
elemType = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

myAssembly.setMeshControls(regions=region, elemShape=TET)
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

# Set unique material properties for each cell

els = myAssembly.instances['pvInstance'].elements
for e in els:

    # Create material based on element's location
    steelTemp = pvModel.Material(name='Steel'+str(e.label))

    # Create the elastic properties: youngsModulus is 209.E3
    # and poissonsRatio is 0.3

    elasticProperties = getElasticProperties(e)
    steelTemp.Elastic(table=(elasticProperties, ) )


    # Create section with region equal to the element region
    
    mySection = pvModel.HomogeneousSolidSection(name='pvSection'+str(e.label),
        material='Steel'+str(e.label), thickness=1.0)

    region = regionToolset.Region(elements=els.sequenceFromLabels([e.label]))
    pvPart.SectionAssignment(region=region,
        sectionName='pvSection'+str(e.label))

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

