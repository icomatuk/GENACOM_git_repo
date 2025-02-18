#!/usr/bin/python
# -*- coding: cp1252 -*-

import sys
import os
from odbAccess import *

job = 'master'
step = 'Step-1'

def fu(job,step):
	
	odb = openOdb(job+'.odb')
	# myNodeSetRF = odb.rootAssembly.instances['PART-1-1'].nodeSets['TOP_REF_POINT']
	# myNodeSetU = odb.rootAssembly.instances['PART-1-1'].nodeSets['BOT_REF_POINT']
	thisStep = odb.steps[step]
	u1 = []
	u2 = []
	u3 = []
	unorm = []
	du1_dt1 = []
	du2_dt1 = []
	du3_dt1 = []
	dunorm_dt1 = []
	
	# Schleife über Frames
	for frameNr in range(0,len(thisStep.frames)):
	    du1_dt1.append( thisStep.frames[frameNr].fieldOutputs['d_U_t_1'].getSubset(position=NODAL).values[59].data[0] )
	    du2_dt1.append( thisStep.frames[frameNr].fieldOutputs['d_U_t_1'].getSubset(position=NODAL).values[59].data[1] )
	    du3_dt1.append( thisStep.frames[frameNr].fieldOutputs['d_U_t_1'].getSubset(position=NODAL).values[59].data[2] )
	    u1.append( thisStep.frames[frameNr].fieldOutputs['U'].getSubset(position=NODAL).values[59].data[0] )
	    u2.append( thisStep.frames[frameNr].fieldOutputs['U'].getSubset(position=NODAL).values[59].data[1] )
	    u3.append( thisStep.frames[frameNr].fieldOutputs['U'].getSubset(position=NODAL).values[59].data[2] )
	    unorm.append( thisStep.frames[frameNr].fieldOutputs['U'].getSubset(position=NODAL).values[59].data[3] )
	    
	odb.close()

	# Ergebnisdatei
	resFile = open(job+'_Fu.txt','w')
	#resFile.write('   u     F\n')
	for i in range(len(unorm)):
	  resFile.write(str(unorm[i])+'\n')
	resFile.close()

def buckling(job,step):

	odb = openOdb(job+'.odb')
	myNodeSetLB = odb.rootAssembly.instances['PART-1-1'].nodeSets['LOCAL_BUCKLE_NODES']
	myNodeSetGB = odb.rootAssembly.instances['PART-1-1'].nodeSets['GLOBAL_BUCKLE_NODE']
	thisStep = odb.steps[step]
	uMaxLB = []
	uGB = []
	
	# Schleife über Frames
	for frameNr in range(0,len(thisStep.frames)):
	  
	  # local buckling
	  thisValues = thisStep.frames[frameNr].fieldOutputs['U'].getSubset(region=myNodeSetLB).values
	  u3=[]
	  for i in range(0, len(thisValues)):
	    u3.append( thisValues[i].data[2] )
	  uMaxLB.append( max(max(u3),-min(u3)) )

	  # global buckling
	  uGB.append( thisStep.frames[frameNr].fieldOutputs['U'].getSubset(region=myNodeSetGB).values[0].data[2] )

	odb.close()

	# Ergebnisdatei
	resFile = open(job+'_GLbuckling.txt','w')
	#resFile.write(' local buckling nodes    global buckling nodes\n')
	for i in range(len(uMaxLB)):
	  resFile.write(str(uMaxLB[i])+'  ,  '+str(uGB[i])+'\n')
	resFile.close()

def materialDamage(job,step):
	
	myFieldOutputKeys = ['HSNFCCRT','HSNFTCRT','HSNMCCRT','HSNMTCRT']
	odb = openOdb(job+'.odb')
	myElementSet = odb.rootAssembly.instances['PART-1-1'].elementSets['SHELLELEM']
	thisStep = odb.steps[step]
	
	damageFrame = [0]*len(myFieldOutputKeys)
	damageElem = [0]*len(myFieldOutputKeys)
	damageIntPoint = [0]*len(myFieldOutputKeys)

	# Schleife über Kriterien
	for fi in range(len(myFieldOutputKeys)):
	  #print 'check '+myFieldOutputKeys[fi]

	  fieldOutputKey = myFieldOutputKeys[fi]

	  # Schleife über Frames
	  for frameNr in range(0,len(thisStep.frames)):
	    #print 'FRAME '+str(frameNr)

	    thisFieldOutput = thisStep.frames[frameNr].fieldOutputs[fieldOutputKey].getSubset(region=myElementSet)

	    # Schleife über alle Values
	    for i in range(0, len(thisFieldOutput.values)):
		# Wenn Kriterium erfüllt ist
		if thisFieldOutput.values[i].data >=1:
		    damageFrame[fi] = frameNr
		    damageElem[fi] = thisFieldOutput.values[i].elementLabel
		    damageIntPoint[fi] = thisFieldOutput.values[i].integrationPoint
		    break # es wird nur ein Wert gesucht, nicht der Maximale

	    if damageFrame[fi]!=0:
		break # wenn der Frame mit dem ersten Überschreiten gefunden ist, brich ab

	odb.close()

	# Ergebnisdatei
	resFile = open(job+'_matDamage.txt','w')
	for fi in range(len(myFieldOutputKeys)):
	  resFile.write(' cirterion: '+myFieldOutputKeys[fi]+'\n')
	  resFile.write('    frame: '+str(damageFrame[fi])+'\n')
	  resFile.write('    element: '+str(damageElem[fi])+'\n')
	  resFile.write('    integration point: '+str(damageIntPoint[fi])+'\n')
	resFile.close()

def cohesiveDamage(job,step):

	fieldOutputKey = 'QUADSCRT'
	odb = openOdb(job+'.odb')
	myElementSet = odb.rootAssembly.instances['PART-1-1'].elementSets['COHESELEM']
	thisStep = odb.steps[step]

	damageFrame = 0
	damageElem = 0
	
	# Schleife über Frames
	for frameNr in range(0,len(thisStep.frames)):
	    #print 'FRAME '+str(frameNr)

	    thisFieldOutput = thisStep.frames[frameNr].fieldOutputs[fieldOutputKey].getSubset(region=myElementSet)

	    # Schleife über alle Values
	    for i in range(0, len(thisFieldOutput.values)):
		# Wenn Kriterium erfüllt ist
		if thisFieldOutput.values[i].data >=1:
		    damageFrame = frameNr
		    damageElem = thisFieldOutput.values[i].elementLabel
		    break # es wird nur ein Wert gesucht, nicht der Maximale

	    if damageFrame!=0:
		break # wenn der Frame mit dem ersten Überschreiten gefunden ist, brich ab

	odb.close()

	# Ergebnisdatei
	resFile = open(job+'_cohDamage.txt','w')
	resFile.write('    frame: '+str(damageFrame)+'\n')
	resFile.write('    element: '+str(damageElem)+'\n')
	resFile.close()

# runtime
fu(job,step)
# buckling(job,step)
# materialDamage(job,step)
# cohesiveDamage(job,step)



# odb = session.odbs[
# odb.steps['Step-1'].frames[50].fieldOutputs['HSNFCCRT'].values[2].data
# odb.steps['Step-1'].frames[50].fieldOutputs['HSNFCCRT'].values[2].elementLabel
# odb.steps['Step-1'].frames[50].fieldOutputs['HSNFCCRT'].values[2].integrationPoint
