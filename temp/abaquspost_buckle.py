# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior


# 后处理
jobname = 'sdy_buckle'

odb = visualization.openOdb(jobname + '.odb')

eigenValueDescription = odb.steps['Step-1'].frames[1].description
eigenValueLoc = eigenValueDescription.index('=') + 1
eigenValue = float(eigenValueDescription[eigenValueLoc:])


outfile = open('buckle_eigenvalue.txt','a')

outfile.write('%f\n' % eigenValue)

outfile.close()
