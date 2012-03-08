"""
Last Modified Aug 10 2011
by Derek Park

This program is designed to analyze growth curves measured by
the 96 well plate TECAN machine. 

The methods are contained within an Analyzer class whose main() method
runs through the process of analyzing the curves.

"""

import os
import sys
import matplotlib
import matplotlib.pyplot as mplot
import numpy
import scipy
from scipy import stats


class Well:

	"""This class encapsualtes all of the relevant information for the well. 
	The alphanumeric label (i.e. 'A1')
	The OD600 measurements for the well
	The strain name in the well
	The dilution in the well"""
	
	
	#Initiliaze a well by passing its label (a string), 
	#the measurements in order (a list), the strain name (a string)
	#the dilution (a string)
	#In addition, there is a variable start_after_drop which is the index
	#of the first measurement after the drop artifact. This is initialized as -1
	#for every Well object, but can be set either from reading from a file 
	#that has the assignments for each drop, or 
	def __init__(self, well_label, msmnts, strain_name, dil):
		self.label = well_label 
		self.measurements = msmnts 
		self.name = strain_name
		self.dilution = dil
		self.starting_timepoint = -1
		self.base = -1
	
	
	def getLabel():
		return self.label
	
	#Return the whole set of measurements
	def getMeasurements(self, start = -1, end = -1):
		if start > -1 and end > -1:
			return self.measurements[start:end]
		else:
			return self.measurements
	
	#Return a specific measurement. i is an integer
	def getSpecificMeasurement(self, i):
		return self.measurements[i]
	

	
	#Return the label of the well
	def getLabel(self):
		return self.label
	
	#Return the name of the strain
	def getStrainName(self):
		return self.name
	
	#Return the dilution
	def getDilution(self):
		return self.dilution
		
	#Return the starting_timepoint, the index where the artifact ends
	def getStartTimepoint(self):
		return self.starting_timepoint
	
	#Set start_after_drop
	def setStartTimepoint(self, i):
		self.starting_timepoint = i	
	
	#Set base
	def setBase(self, base):
		self.base = float(base)
	
	#Get base
	def getBase(self):
		return self.base
		
	#Get the measurements minus the base
	def getMeasurementsLessBase(self):
			
		baseLessMeasurements = []
		
		n = 0
		while n < len(self.measurements):
			baseLessMeasurements.append(float(self.measurements[n]) - self.base)
			n += 1
		return baseLessMeasurements

class Analyzer:
	
	"""The Analyzer class loads the data and has methods for analyzing it."""
	
	def __init__(self, OD_FILE_NAME, LABEL_FILE_NAME):
	
		self.OD600_WELLS, self.timesHrs, self.timesSecs = self.load_OD600(OD_FILE_NAME, LABEL_FILE_NAME)
		self.doublingTimes = {}
		
		
	####################################################################################################################
	#   This method loads the OD600 values and the well labels. It then constructs individual Well classes for each well
	#   Since this data is read from file, the formatting needs to be correct. 
	#
	#
	#   For the data file that has the OD600 file, it should be TAB DELIMITED. It should also be a copy of the OD600
	#   Section of the TECAN data file. For example:
	#
	#	OD600
	#	Cycle Nr.		1		2		3		.....
	#	Time [s]		0		300		600		.....
	#	Temp [C]		37		37		37		.....
	#	A1				0.01	0.02	0.03	.....
	#	A2				0.05	0.01	0.10	.....
	#	.
	#	.
	#	.
	#	H12				0.01	0.02	0.03	.....
	#
	#
	#
	#   For the data file that has the labels for the different wells, it needs to be hand made, but should be formatted 
	#   as such:
	#
	#	A1		[Strain Name]		[Dilution]
	#	A2		[Strain Name]		[Dilution]
	#	A3		[Strain Name]		[Dilution]
	#	.
	#	.
	#	.
	#	H12 	[Strain Name]		{Dilution]
	#
	#
	#   The method returns 3 lists:
	#   Wells, a list of the wells with their corresponding measurements/data
	#   Hours, a list of the measurement times converted into hours
	#   Seconds a list of the measurement times in seconds
	####################################################################################################################
	
	def load_OD600(self, OD600s, labels):
		Wells = [] #a list of Wells
		Hours = [] #A list of floats
		Seconds = [] #A list of floats
		
		
		#Open the file of OD600 data
		fileOD600 = open(OD600s)
		
		fileOD600.readline() #Skip the first line, it's just the 'OD600' label
		fileOD600.readline() #Skip the second line, it's just the cycle numbers
		
		#Read the next line which has the timepoints in seconds. Add that to Seconds and also convert the times
		#into hours and add them to Hours
		
		timePointsInSeconds = fileOD600.readline().split("\t")
		
		n = 1 #We start at index = 1 (i.e. the second) token in this line because the first one is just the string label 'Times [s]'
		while n < len(timePointsInSeconds):
			currTimePointSeconds = float(timePointsInSeconds[n])
			currTimePointHours = float(currTimePointSeconds/3600)
			
			Hours.append(currTimePointHours)
			Seconds.append(currTimePointSeconds)
			
			n+=1
		
		fileOD600.readline() #Skip the next line because it is just the temperatures
		
		#Open the file of well label annotations. NOTE: THE WELLS SHOULD BE SORTED FROM A1 TO H12. THAT WAY IT MATCHES THE OD600 DATA ORDER
		fileLabels = open(labels)
		
		
		#Read each well and create a new Well for each one and add it to the Wells list
		n = 0
		while n < 96:
			currOD600Line = fileOD600.readline() #Read each line
			currLabelLine = fileLabels.readline() 
			
			currOD600Tokens = currOD600Line.split("\t")
			currLabelTokens = currLabelLine.split("\t")
			
			currLabel = currOD600Tokens[0]
			currMeasurements = currOD600Tokens[1:len(currOD600Tokens)]
			currStrainName = currLabelTokens[1]
			currDilution = currLabelTokens[2]
			
			currWell = Well(currLabel, currMeasurements, currStrainName, currDilution)
			Wells.append(currWell)
			
			n +=1
			
		
		#close the files
		fileOD600.close()
		fileLabels.close()
		
		
		return Wells, Hours, Seconds
	
	
	
	
	#Plots the data for a given well name. Parameters specifying the starting index of the data and the ending index of the data can be
	#specified, but are optional. The default is that the whole series is plotted
	#X axis will be hours
	def plotWell(self, well_name, start = -1, end = -1):
		
		#First, search through the wells and find the specified well 
		currWell = None
		n = 0
		while n < len(self.OD600_WELLS):
			
			possibleWell = self.OD600_WELLS[n]
			 
			if possibleWell.getLabel() == well_name:
				currWell = possibleWell
				n = len(self.OD600_WELLS)
			else:
				n +=1
		
		#now get the measurements to be plotted
		dataToBePlotted = currWell.getMeasurements(start,end) #The exception case where start and end are not specified is handled by the getMeasurements method
		
		
		#now get the x axis values corresponding to the measurements
		xAxisValues = []
		if start == -1 and end == -1:
			xAxisValues = self.timesHrs
		else:
			#a start and an ending index have been specified
			xAxisValues = self.timesHrs[start:end] 
		
		
		#Finally, plot it
		mplot.plot(xAxisValues,dataToBePlotted)
	
	
	
	# This method tries to load the starting timepoints from the file name specified. 
	# The file must be in this format:
	#
	#
	# A1		[Index]
	# A2		[Index]
	# .
	# .
	# .
	# H12		[Index]
	#
	#
	# File should be tab delimited. 
	# In the event that no such file exists, it will read the corresponding name and index,
	# find the corresponding well and define it with the index.
	# 
	# If there is no file, it calls createStartTimepoints() to build the list.
	def loadStartTimepoints(self,fileName):
		
		try:
			file = open(fileName)
			
			for line in file:
				currLineEntries = line.split("\t")
				currWellName = currLineEntries[0]
				currIndex = int(currLineEntries[1])
				
				#Now loop through the list of wells, find the right well, and assign its index
				n = 0
				while n < len(self.OD600_WELLS):
					currWell = self.OD600_WELLS[n]
					
					if currWell.getLabel() == currWellName:
						currWell.setStartTimepoint(currIndex)
						n = len(self.OD600_WELLS)
					
					n += 1
				
			print "Indeces read and loaded from file ", fileName	
				
		except IOError:
			print "Failure to load file. Creating new one"
			self.createStartTimepoints(fileName)
	
	
	
	# This method plots the beginning of each growth curve for each well
	# and gives a guess as to where the drop is.
	# 
	# It calculates the max in the range and the timepoint at which it occurs.
	# If it is before the end, i.e. the max occurs in the middle of the window, 
	# Then it suggests that it is an artifact. If it is an artifact, the index is 
	# set to the idnex after the peak. Else, it is 0, i.e. the beginning. 
	# 
	# Alternatively the user can manually set the post-drop index. They can review the 
	# data which is printed to the console
	#
	# User input controls the decision. At the end, the data is written to the file 'fileName'
	def createStartTimepoints(self, fileName):
		
		#These determine the window of values to be plotted. 0 means beginning.
		startTimepoint = 0
		endTimepoint = 100
		
		
		#First, loop through the wells and plot each well, asking for user input
		#about how to set the beginning timepoint
		n = 0
		while n < len(self.OD600_WELLS):
			
			currWell = self.OD600_WELLS[n]
			print "N is: ", n, "  !!!!!!!!!!!!"
			#get the subrange of X and Y values to be plotted
			currODVals = currWell.getMeasurements(startTimepoint, endTimepoint)
			currTimepoints = range(startTimepoint, endTimepoint)
		
			#get the maximum OD in this range as well as the timepoint at which it 
			#occurs
			maxOD = max(currODVals)
			maxTimepoint = currODVals.index(maxOD)
			
			#Print out the range of values, adding a marker to show the max
			i = 0
			print "\n\nData for well ", currWell.getLabel()
			print "Timepoint		OD600"
			print "----------------------"
			
			while i < len(currODVals):
				
				if i == maxTimepoint:
					print currTimepoints[i], "   ", currODVals[i], " <----- MAX"
				else:
					print currTimepoints[i], "   ", currODVals[i]
				
				i +=1
			
			
			
			
			#Plot the values
			mplot.plot(currTimepoints, currODVals, "b+")
			mplot.title("Well " + currWell.getLabel() + " " + currWell.getStrainName() + " " + currWell.getDilution())
			mplot.xlabel("Timepoints")
			mplot.ylabel("OD600")
			
			
			
			#Ask for user input until a valid input is given
			seekInput = 1
			while seekInput:
				print "\n\n\nData for well ", currWell.getLabel()
				print "Current maximum is: ", maxOD,
				print "This occurs at timepoint: ", maxTimepoint
				print "Options:"
				print "'y'      - Accept this timepoint as beginning timepoint"
				print "'n'      - Manually enter beginnign timepoint"
				print "'resize' - resize the graph window" 
				print "'exit'   - Exit this method"
				
				input = raw_input()
				
				if input == 'exit':
					seekInput = 0
					n = len(self.OD600_WELLS)
				elif input == 'y':
					seekInput = 0
					currWell.setStartTimepoint(maxTimepoint +1)
					print "Beginning timepoint for well ", currWell.getLabel(), " set to ", (maxTimepoint+1)
				elif input == 'n':
					print "Input new starting timepoint"
					newTimepoint = int(raw_input())
					currWell.setStartTimepoint(newTimepoint)
					print "Beginning timepoint for well ", currWell.getLabel(), " set to ", newTimepoint
					seekInput = 0
				elif input == "resize":
					print "Currently displaying timepoints from timepoint ", startTimepoint, " to timepoint ", endTimepoint
					print "Enter new starting timepoint"
					startTimepoint = int(raw_input())
					print "Enter new ending timepoint"
					endTimepoint = int(raw_input())
					n -= 1
					seekInput = 0
			
			n +=1
			mplot.close()
			
			
		
		
		n = 0
		
		
		while n < len(self.OD600_WELLS):
			currWell = self.OD600_WELLS[n]
			#Now save it to the file by appending it
			f = open(fileName, 'a')
			
			f.write(currWell.getLabel() + "\t")
			f.write(str(currWell.getStartTimepoint()))
			f.write("\n")
			f.close()
			n +=1
			
	
	
	
	# This method finds the base for each well.
	# The base is defined as the minimum of the first 10 values
	# From the starting timepoint. The starting timepoint is 
	# got by calling the Well object's getStartTimepoint() method.
	def findBase(self):
		
		n = 0
		while n < len(self.OD600_WELLS):
			currWell = self.OD600_WELLS[n]
			
			currStartTimepoint = currWell.getStartTimepoint()
			currEndTimepoint = currStartTimepoint + 10 
			
			currWellVals = currWell.getMeasurements(currStartTimepoint, currEndTimepoint)

			currWell.setBase(min(currWellVals))
			print currWell.getLabel(), " Starting index is: ",currStartTimepoint, " base is: ", currWell.getBase()
			n +=1 
		
	

	# This method finds all of the doubling times for each well using a specified windowSize
	# The results will be returned as a dictionary with the key being the alphanumeric
	# label of the well and the entry being a list of doubling times. 
	# The list will have all of the doubling times from timepoint 0.
	# This includes artifactual timepoints with the spike at the beginning involved which can later be excluded
	def findDoublingTimes(self, windowSize):
		doublings = {}
		
		wellIndex = 0
		while wellIndex < len(self.OD600_WELLS):
			
			currDoublings = []
			
			currWell = self.OD600_WELLS[wellIndex]
			
			
			#Get the measurements for this well MINUS THE BASE
			currMeasurements = currWell.getMeasurementsLessBase()
			
			#The problem with currMeasurements, though, is that it may contain negative values.
			#This poses a problem when we want to log transform the data because a log transform of
			#a number <= 0 results in either NaN or -inf, both of which will make our linear regressions
			#problematic (infinite loop). 
			#
			#The solution is to go through the list, find the 0 or negative values, and replace them
			#If such a value is at timepoint 0, we replace it with 0.000000001, or 10E-9, a close enough
			#approximation of 0. If the value occurs in the middle of the list, (timepoint > 0), then the
			#the measurement at timepoint n is equal to the measurement at timepoint n-1.
			adjustedMeasurements = []
			
			
			#First, test if the first measurement is <= 0, if so, replace it with 10E-9
			firstMeasurement = currMeasurements[0]
			if firstMeasurement <= 0:
				adjustedMeasurements.append(0.000000001)
			else:
				adjustedMeasurements.append(firstMeasurement)
			
			#Now go through the rest of the list
			measurementIndex = 1
			while measurementIndex < len(currMeasurements):
				
				toBeAdded = 0
				
				possibleAddition = currMeasurements[measurementIndex]
				
				if possibleAddition <= 0:
					toBeAdded = adjustedMeasurements[measurementIndex-1]
				else:
					toBeAdded = possibleAddition
				
				adjustedMeasurements.append(toBeAdded)
				
				
				measurementIndex +=1
			
			
			#Now log transform (base 2) the measurements:
			logMeasurements = numpy.log2(adjustedMeasurements)
			
			#Now we calculate the doubling rate using a sliding window
			startOfWindow = 0
			while (startOfWindow + windowSize) < len(logMeasurements):
				
				currYVals = logMeasurements[startOfWindow:(startOfWindow + windowSize)]
				currXVals = self.timesHrs[startOfWindow: (startOfWindow + windowSize)]
				
				currSlope,currIntercept, currR, currP, currErr = scipy.stats.linregress(currXVals, currYVals)
				
				#print currSlope
				currDoublingTime = 0
				if currSlope != 0:
					currDoublingTime = 1/float(currSlope)
				
				currDoublings.append(currDoublingTime)
				
				
				startOfWindow +=1
			
			
			
			#Now that we have the doublings associated with this well, add it to the dictionary
			currLabel = currWell.getLabel()
			
			doublings[currLabel] = currDoublings
			
			
			wellIndex += 1
		
		self.doublingTimes = doublings
		return doublings
	

	#Saves to file with name fileName
	#The file will be tab delimited and have this format:
	# OD600 Doubling Times
	# -			-			-				-										Timepoint:		0					1				...
	# Well		Strain		Dilution		First Timpepoint After Artifact			Interval: 		[interval 1]		[interval 2]	...\
	# A1
	# A2
	# .
	# .
	# .
	def saveToFile(self, fileName, windowSize):
		f = open(fileName, "w")
		f.write("OD600 Doubling Times \n")
		
		timepointLine = "- \t - \t - \t - \t Timepoint: \t"
		n = 0
		while (n+windowSize) < len(self.timesHrs):
			timepointLine += str(n) + "\t"
			n +=1
		
		timepointLine += "\n"
		
		f.write(timepointLine)
		
		thirdLine = "Well \tStrain \tDilution \tFirst Timepoint After Artifact \tInterval (hrs): \t"
		
		n = 0
		while (n+windowSize) < len(self.timesHrs):
			
			currIntervalStart = str(self.timesHrs[n])
			currIntervalEnd = str(self.timesHrs[(n + windowSize)])
			
			thirdLine += currIntervalStart
			thirdLine += " - "
			thirdLine += currIntervalEnd
			thirdLine += "\t"
			
			n +=1
		
		thirdLine += "\n"
			
		f.write(thirdLine)
		
		
		wellIndex = 0
		while wellIndex < len(self.OD600_WELLS):
			currWell = self.OD600_WELLS[wellIndex]
			
			line = ""
			
			line += str(currWell.getLabel())
			line += "\t"
			
			line += str(currWell.getStrainName())
			line += "\t"
			
			line += str(float(currWell.getDilution()))
			line += "\t"
			
			line += str(currWell.getStartTimepoint())
			line += "\t"
			line += "\t"
			
			currDoublings = self.doublingTimes[currWell.getLabel()]
			
			doublingListIndex = 0
			while doublingListIndex < len(currDoublings):
				line += str(currDoublings[doublingListIndex])
				line += "\t"
				doublingListIndex +=1
			
			line += "\n"
			
			f.write(line)
			
			wellIndex +=1
		
		f.close()
			
			
			

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#Start of the main script
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
def run():

	"""Global Parameters
		
	Define the global parameters here. Change these parameters (or, if you're a real do-gooder, 
	make a user-input based definition system) to control which files the Data are read from
	"""
		
		
	#This is the COMPLETE location of the tab delimited text file for the OD600 data values.
	#For information on proper formatting of this file, see the documentation for the load_OD600() method
	OD_600_DATA_FILE_NAME = "/Users/yyfwuhan/Projects/2011-TaMaRa-growth-curve/TECANWellAnalyzer/EXAMPLE/OD600_Values.txt"
		
	#This is the location of the tab delimited text file for the strain name and concentration for the wells
	#See the annotation in the load_OD600() method for how to properly format the file
	WELL_LABEL_FILE_NAME = "/Users/yyfwuhan/Projects/2011-TaMaRa-growth-curve/TECANWellAnalyzer/EXAMPLE/well-labels_7-14-11.txt"

	#This is the location of the file that has the measurement indeces at which the artificial peak at the beginnign
	#of growth curves ends
	DROP_INDECES_FILE_NAME = "/Users/yyfwuhan/Projects/2011-TaMaRa-growth-curve/TECANWellAnalyzer/EXAMPLE/Starts_After_Drops.txt"


	#This is the window size when calculating the doubling time. 
	DOUBLING_WINDOW_SIZE = 40

	#This is the file name for printing the doubling times to file
	DOUBLING_FILE_NAME = "/Users/yyfwuhan/Projects/2011-TaMaRa-growth-curve/TECANWellAnalyzer/EXAMPLE/Doubling_time.txt"

	a = Analyzer(OD_600_DATA_FILE_NAME, WELL_LABEL_FILE_NAME)

	#First: Try to load the indeces ater drop
	a.loadStartTimepoints(DROP_INDECES_FILE_NAME)

	#Second: Find the base OD for each well
	a.findBase()

	#Third: Calculate the complete set of doubling times for each well. It will be returned as a dictionary
	doublingTimes = a.findDoublingTimes(DOUBLING_WINDOW_SIZE)	

	a.saveToFile(DOUBLING_FILE_NAME, DOUBLING_WINDOW_SIZE)

