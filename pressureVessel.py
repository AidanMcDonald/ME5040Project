"""
pressureVessel.py

Model a simple reactor pressure vessel with continuously
changing material properties.
"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import part
import inspect
import material
import section
import assembly
import step
import load
import mesh
import job
import visualization
from math import log


backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)


def printMethods(object):
    # Abaqus has terrible documentation so this function is useful for debugging
    object_methods = [method_name for method_name in dir(object)]

def getMaterialProperties(element, part, sourceLocation, t):
    # Set material properties for an element
    # t in years
    
    # Get location of element from nodes
    elementNodeLabels = element.connectivity
    nodes = part.nodes
    coords = []
    for n in elementNodeLabels:
        coords.append(nodes.getFromLabel(n+1).coordinates)
    centroid = [0.0, 0.0, 0.0]
    for i in range(3):
        centroid[i] = sum([c[i] for c in coords])/len(elementNodeLabels)
    
    # Estimate neutron radiation from location
    distanceFromSource = sum([(centroid[i]-sourceLocation[i])**2 for i in range(3)])**.5  # [m]
    Phi_reactor = 1.25e18  # [1/(cm^2*yr)]
    Phi = Phi_reactor*1.52**2/distanceFromSource**2  # [1/(cm^2*yr)]
    
    sigmaBar_d = 249.24*1e-24  # [cm^2]
    dpa = Phi*t*sigmaBar_d  # []
    if dpa>3.0:
        dpa=3.0
    
    # Set material properties based on dpa
    youngsModulus = 209E9
    poissonsRatio = 0.3
    if dpa==0.0:
        sigma_y = 300e6  # [Pa]
    else:
        sigma_y = (52.93*log(dpa)+729.5)*1e6  # [Pa]
    
    epsilon_y = sigma_y/youngsModulus  # []
    
    # Interpolate the values from the graph in the paper to get sigma_ut, epsilon_ut

    dpa_table = [0.0, 0.4, 1.25, 3.0]
    UT_table = [500.0e6, 680.0e6, 760.0e6, 790.0e6]
    UE_table = [.48, .1, .06, .005]
    
    for i in range(len(dpa_table)):
        if dpa>=dpa_table[i] and dpa<=dpa_table[i+1]:
            interpBin=i
            break
    
    x_interp = (dpa-dpa_table[interpBin])/(dpa_table[interpBin+1]-dpa_table[interpBin])
    sigma_ut = UT_table[interpBin]+x_interp*(UT_table[interpBin+1]-UT_table[interpBin])
    epsilon_ut = UE_table[interpBin]+x_interp*(UE_table[interpBin+1]-UE_table[interpBin])  
    
    elasticProperties = (youngsModulus, poissonsRatio)
    plasticProperties = ((sigma_y, epsilon_y-epsilon_y),(sigma_ut, epsilon_ut-epsilon_y))
    
    return elasticProperties, plasticProperties


def ExecutePressureVesselCase(modelName, t, pressure, sourceLocation, meshSize):
    # Get model constants for AP1000
    innerRadius = 4.0386/2  # [m]
    thickness = .203  # [m]
    designPressure = 17.2e6  # [Pa(a)]
    designTemperature = 343.3  # [C]
    height = 12.056  # [m]
    
    # Create a model.
    pvModel = mdb.Model(name=modelName)

    # Create a new viewport in which to display the model
    # and the results of the analysis.

    myViewport = session.Viewport(name='Pressure Vessel',
        origin=(20, 20), width=150, height=120)
        
    #-----------------------------------------------------
    # Create a sketch for the base feature.

    pvRotateSketch = pvModel.ConstrainedSketch(name='pvRotateSketch',
        sheetSize=100.)

    pvRotateSketch.ArcByCenterEnds(center=(0,height/2-innerRadius), 
                                   point2=(0, height/2),
                                   point1=(innerRadius, height/2-innerRadius))
    pvRotateSketch.Line(point1=(innerRadius, height/2-innerRadius),
                        point2=(innerRadius, -height/2+innerRadius))
    pvRotateSketch.ArcByCenterEnds(center=(0, -height/2+innerRadius),
                                   point2=(innerRadius, -height/2+innerRadius),
                                   point1=(0, -height/2))
    pvRotateSketch.Line(point1=(0, -height/2),
                        point2=(0, -height/2-thickness))
    pvRotateSketch.ArcByCenterEnds(center=(0, -height/2+innerRadius),
                                   point1=(0, -height/2-thickness),
                                   point2=(innerRadius+thickness, -height/2+innerRadius))
    pvRotateSketch.Line(point1=(innerRadius+thickness, -height/2+innerRadius),
                        point2=(innerRadius+thickness, height/2-innerRadius))
    pvRotateSketch.ArcByCenterEnds(center=(0,height/2-innerRadius),
                                   point1=(innerRadius+thickness, height/2-innerRadius),
                                   point2=(0, height/2+thickness))
    pvRotateSketch.Line(point1=(0, height/2+thickness),
                        point2=(0, height/2))
    pvRotateSketch.assignCenterline(pvRotateSketch.ConstructionLine(point1=(0,0), point2=(0,1)))

    sketch = mdb.models[modelName].convertAllSketches()

    # Create a three-dimensional, deformable part.

    pvPart = pvModel.Part(name='PressureVessel', dimensionality=THREE_D,
        type=DEFORMABLE_BODY)

    # Create the part's base feature by extruding the sketch 
    # through a distance of 25.0.

    pvPart.BaseSolidRevolve(sketch=pvRotateSketch, angle=360)


    #-----------------------------------------------------
    # Create a part instance.

    myAssembly = pvModel.rootAssembly
    myInstance = myAssembly.Instance(name='pvInstance',
        part=pvPart, dependent=ON)

    #-------------------------------------------------------
    # Create a step. The time period of the static step is 1.0, 
    # and the initial incrementation is 0.1; the step is created
    # after the initial step. 

    pvModel.StaticStep(name='pressureStep', previous='Initial',
        timePeriod=1.0, initialInc=0.1,
        description='Pressurize it.')

    #-------------------------------------------------------
    # Create a pressure load.
    insideFacePoint = (innerRadius,0,0)
    insideFace = myInstance.faces.findAt((insideFacePoint,) )

    insideTopFacePoint = (0,height/2,0)
    insideTopFace = myInstance.faces.findAt((insideTopFacePoint,) )

    insideBottomFacePoint = (0,-height/2,0)
    insideBottomFace = myInstance.faces.findAt((insideBottomFacePoint,) )

    insideSurfaces = ((insideFace, SIDE1), (insideTopFace, SIDE1), (insideBottomFace, SIDE1))

    pvModel.Pressure(name='Pressure', createStepName='pressureStep',
        region=insideSurfaces, magnitude=pressure)

    #-------------------------------------------------------
    # Create a boundary condition so it doesn't move or rotate
    bottomPoint = myInstance.vertices.findAt(((0,-height/2-thickness,0),))
    pvModel.DisplacementBC(name='PinAtBottom', createStepName = 'Initial', region=(bottomPoint,),
                           u1=0,u2=0,u3=0,ur1=0,ur2=0,ur3=0)
    topPoint = myInstance.vertices.findAt(((0,height/2+thickness,0),))
    pvModel.DisplacementBC(name='PinAtTop', createStepName = 'Initial', region=(topPoint,),
                           u1=0,u3=0,ur1=0,ur2=0,ur3=0)

    #-------------------------------------------------------
    # Assign an element type to the part instance.

    region = [cell for cell in pvPart.cells]
    elemType = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

    pvPart.setMeshControls(regions=region, elemShape=TET)
    pvPart.setElementType(regions=region, elemTypes=(elemType,))

    # Seed the part instance.

    pvPart.seedPart(size=meshSize)

    # Mesh the part instance.

    pvPart.generateMesh()

    # Display the meshed beam.

    myViewport.assemblyDisplay.setValues(mesh=ON)
    myViewport.assemblyDisplay.meshOptions.setValues(meshTechnique=ON)
    myViewport.setValues(displayedObject=pvPart)

    #-------------------------------------------------------
    # Set unique material properties for each cell
    els = pvPart.elements
    for e in els:

        # Create material based on element's location
        steelTemp = pvModel.Material(name='Steel'+str(e.label))

        elasticProperties, plasticProperties = getMaterialProperties(e, pvPart, sourceLocation, t)
        steelTemp.Elastic(table=(elasticProperties, ) )
        steelTemp.Plastic(table=plasticProperties)

        # Create section with region equal to the element region
        mySection = pvModel.HomogeneousSolidSection(name='pvSection'+str(e.label),
            material='Steel'+str(e.label), thickness=.25)

        region = regionToolset.Region(elements=els.sequenceFromLabels([e.label]))
        pvPart.SectionAssignment(region=region,
            sectionName='pvSection'+str(e.label))

    #-------------------------------------------------------
    # Create an analysis job for the model and submit it.

    jobName = modelName
    myJob = mdb.Job(name=jobName, model=modelName,
        description='pressure vessel project')

    # Wait for the job to complete.

    myJob.submit()
    myJob.waitForCompletion()

    #-------------------------------------------------------
    # Open the output database and display a
    # default contour plot.

    myOdb = visualization.openOdb(path=jobName + '.odb')
    myViewport.setValues(displayedObject=myOdb)
    
if __name__ == '__main__':
    innerRadius = 4.0386/2  # [m]
    thickness = .203  # [m]
    designPressure = 17.2e6  # [Pa(a)]
    designTemperature = 343.3  # [C]
    height = 12.056  # [m]
    
    #ExecutePressureVesselCase(modelName="ControlCase", t=0, pressure=1.5e7, sourceLocation = [0.0, (-height/2+innerRadius)/2, 0.0], meshSize=1.0)
    #ExecutePressureVesselCase(modelName="NominalCase", t=60, pressure=1.5e7, sourceLocation = [0.0, (-height/2+innerRadius)/2, 0.0], meshSize=1.0)
    #ExecutePressureVesselCase(modelName="ReallyCookIt", t=6e4, pressure=1.5e7*5, sourceLocation = [0.0, (-height/2+innerRadius)/2, 0.0], meshSize=1.0)
    #ExecutePressureVesselCase(modelName="OutsideSource", t=6e4, pressure=1.5e9, sourceLocation = [innerRadius+thickness+.2, (-height/2+innerRadius)/2, 0.0], meshSize=1.0)

    ExecutePressureVesselCase(modelName="DontPopIt", t=0, pressure=1.5e7*5.5, sourceLocation = [0.0, (-height/2+innerRadius)/2, 0.0], meshSize=1.0)
    ExecutePressureVesselCase(modelName="PopIt", t=6e4, pressure=1.5e7*5.5, sourceLocation = [innerRadius+thickness+.2, (-height/2+innerRadius), 0.0], meshSize=1.0)