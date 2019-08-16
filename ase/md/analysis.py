from scipy import stats
from math import floor
from random import randint
from collections import Counter
import numpy as np
import builtins
import os

import matplotlib.pyplot as plt 

from ase.io import read, write, ulm
from ase.io.trajectory import Trajectory
from ase import units


def frange(start, stop, step=1.0):
    ''' "range()" like function which accept float type''' 
    i = start
    while i < stop:
        yield i
        i += step

class DiffusionCoefficient:
	def __init__(self, file_name, timestep, steps_between_saved_images, no_of_segments=1, ignore_n_images=0, continuous=True, molecule_index=None, plot=False, data=False):
		'''Caclulates Diffusion Coefficients and Plots the graph used for calculation

		Parameters:
			file_name (String): Name of the .traj file where the images are stored
			timestep (Int): Timestep used between each images
			steps_between_saved_images (Int): Interval used when writing the .traj file 
			no_of_segments (Int): Divides the given images to segments to prevent from loading all the images in the RAM at once
			ignore_n_images (Int): Number of images you want to ignore from the start for which you don't want to calculate
			continuous (Bool): Set True if the user wants the graph to be continuous
			molecule_index (List of Int): The indices of atoms whose Diffusion Coefficient is to be calculated
			plot (Bool): If the user wants to plot the graph or not
			data (Bool): If the user wants to print the data or not

			This function calculates the Diffusion Coefficient for the given .traj file using the Einstein Equation:
			⟨|r(t)−r(0)|**2⟩=6Dt (where r(t) is the position of atom at time t, D is the Diffusion Coefficient)
			wiki : https://en.wikibooks.org/wiki/Molecular_Simulation/Diffusion_Coefficients
		'''

		self.input_file = file_name
		self.timestep = timestep
		self.steps_between_saved_images = steps_between_saved_images
		self.no_of_segments = no_of_segments
		self.ignore_n_images = ignore_n_images
		self.continuous = continuous
		self.molecule_index = molecule_index
		self.plot_graph = plot
		self.data = data

		# Getting information for types of atoms from given file so as to calculate diffusion coefficients of different types of atoms separately
		temp=read(filename=self.input_file, index='0:1')
		write(filename='temp.traj',images=temp)
		traj=Trajectory('temp.traj')
		symbols=[]
		for i in range(len(traj[0])):
		    symbols.append(traj[0].symbols[i])
		
		# Condition used if user wants to calculate Diffusion Coefficients for a specific atom or group of atoms (molecule)
		if self.molecule_index == None:
			self.types_of_atoms = list(Counter(symbols).keys())
			self.no_of_types_of_atoms = len(self.types_of_atoms)
			self.no_of_atoms = list(Counter(symbols).values())
		else:
			self.no_of_types_of_atoms = 1
			self.no_of_atoms = 1
			name = ''
			for i in self.molecule_index:
				name += symbols[i]
			self.types_of_atoms = [name]
		
		self.segment_ensemble_average = []
		self.xyz_segment_ensemble_average = []
		for i in range(self.no_of_segments):
			self.xyz_segment_ensemble_average.append([])
			self.segment_ensemble_average.append([])
			for j in range(self.no_of_types_of_atoms):
				self.xyz_segment_ensemble_average[i].append([])
				self.segment_ensemble_average[i].append([])
				for k in range(3):
					self.xyz_segment_ensemble_average[i][j].append([])
		
		self.time_between_images = self.timestep * self.steps_between_saved_images
		
		
 		# To find out total number of images in .traj file
		self.total_images = ulm.read_header(builtins.open(self.input_file,'rb'))[2]
		self.total_images = self.total_images - self.ignore_n_images

		# List of Timestep used for plotting (x-axis)
		self.timesteps = list(frange(0,self.total_images*self.time_between_images,self.time_between_images))

		self.slopes=[]
		self.intercepts=[]
		for i in range(self.no_of_types_of_atoms):
			self.slopes.append([])
			self.intercepts.append([])

		self.cont_xyz_segment_ensemble_average = 0
		self.segment_line=[]
		self.XYZ=['X','Y','Z']
		self.xyz_markers = {'X':'o', 'Y':'s','Z':'^'}
		self.color_list = plt.cm.Set3(np.linspace(0, 1, self.no_of_types_of_atoms))

	def calculate(self):
		for segment_no in range(self.no_of_segments):
			start = self.ignore_n_images + (segment_no*floor(self.total_images/self.no_of_segments))  # Making separate .traj files for each segment, so as to prevent from loading all the images in the RAM at once
			end = start + floor(self.total_images/self.no_of_segments)
			temp = read(filename=self.input_file,index='%d:%d'%(start,end))
			write(filename='temp.traj', images=temp)
			traj = Trajectory('temp.traj') 
		
		
			for image_no in range(1,len(traj)): ##1,len(traj)
				#xyz_ensemble_average = [[0.0]*3]*no_of_types_of_atoms
				xyz_ensemble_average = []
				for i in range(self.no_of_types_of_atoms):
					xyz_ensemble_average.append([0.,0.,0.])
				
				
				if self.molecule_index == None: # Calculating for all the atoms
					for atom_no in range(len(traj[image_no])):
						index = self.types_of_atoms.index(traj[image_no].symbols[atom_no])
						for i in range(3):
							xyz_ensemble_average[index][i] += np.square(traj[image_no].positions[atom_no][i] - traj[0].positions[atom_no][i])
					for xyz in range(3):
						for index in range(self.no_of_types_of_atoms):
							self.xyz_segment_ensemble_average[segment_no][index][xyz].append(xyz_ensemble_average[index][xyz]/(self.no_of_atoms[index]*2))
					for index in range(self.no_of_types_of_atoms):
						self.segment_ensemble_average[segment_no][index].append(np.sum(xyz_ensemble_average[index])/(6*self.no_of_atoms[index]))
		
				else: # Calculating for specific atom or group of atoms (molecule)
					com_orig = np.zeros(3)
					com_new = np.zeros(3)
					for atom_no in self.molecule_index:
						for xyz in range(3):
							com_orig[xyz] += traj[0].positions[atom_no][xyz]
							com_new[xyz] += traj[image_no].positions[atom_no][xyz]
					com_orig /= len(self.molecule_index)
					com_new /= len(self.molecule_index)
					for xyz in range(3):
						xyz_ensemble_average[0][xyz] += np.square(com_new[xyz] - com_orig[xyz])
						self.xyz_segment_ensemble_average[segment_no][0][xyz].append(xyz_ensemble_average[0][xyz]/2)
					self.segment_ensemble_average[segment_no][0].append(np.sum(xyz_ensemble_average[0][xyz]/6))

		self.cont_xyz_segment_ensemble_average = self.xyz_segment_ensemble_average
		os.remove('temp.traj')

	def fit_data_and_plot(self,x, y, plt, color, label=None, cont=False, cont_y=None):
		# Generated linear fit  
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		line = slope*np.asarray(x)+intercept
		if label[0] != '_':
			skip = label.index('A')
		else:
			skip = 0
		mark = self.xyz_markers[label[-1]]

		if cont == False:
			plt.plot(x, y, color=color, marker=mark,label=label[skip:], linewidth=0)
			intercept_c=intercept
		else:
			slope, intercept_c, r_value, p_value, std_err = stats.linregress(x,cont_y)
			line_c = slope*np.asarray(x)+intercept_c
			plt.plot(x, cont_y, color=color, marker=mark,label=label[skip:], linewidth=0)

		if label[0] == '_':
			label=label[1:]

		# print(r'Intercept = %.10f $\AA^2$; Diffusion Coefficient = %.10f $\AA^2$/fs' % (intercept, slope))
		# Converting to more common units of cm^2/s => multiply by 10^-1
		# \AA^2 => cm^2 requires multiplying by (10^-8)^-2
		# fs => s requires dividing by 10^-15
		if self.data == True:
			print('---')
			print(r'%10s: Intercept = %.10f cm^2; Diffusion Coefficient = %.10f cm^2/s, %.10f m^2/s' % (label, intercept_c/(10**8), slope*(0.1), slope*(10**-5)))
		return slope, intercept_c

	def plot(self):
		self.calculate()
		for segment_no in range(self.no_of_segments):
			if segment_no>0 and self.continuous==True:
			    for index in range(self.no_of_types_of_atoms):
			        for xyz in range(3):
			            self.cont_xyz_segment_ensemble_average[segment_no][index][xyz]=np.array(self.cont_xyz_segment_ensemble_average[segment_no][index][xyz])+self.cont_xyz_segment_ensemble_average[segment_no-1][index][xyz][-1]
			    
			for index in range(self.no_of_types_of_atoms):	
				for xyz in range(3):
					if segment_no == 0:
						custom_label = 'Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], self.XYZ[xyz])
					else:
						custom_label = '_Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], self.XYZ[xyz]) # To remove label for other segments as all the segments are of same colour
					ratio = floor(self.total_images/self.no_of_segments)
					slope, intercept = 	self.fit_data_and_plot(self.timesteps[segment_no*ratio+1:segment_no*ratio+ratio], self.xyz_segment_ensemble_average[segment_no][index][xyz], plt, self.color_list[index], custom_label, self.continuous, self.cont_xyz_segment_ensemble_average[segment_no][index][xyz])
					self.slopes[index].append(slope)
					self.intercepts[index].append(intercept)

			self.segment_line.append(self.timesteps[segment_no*floor(self.total_images/self.no_of_segments)+floor(self.total_images/self.no_of_segments)-1] + self.time_between_images)


		for index in range(self.no_of_types_of_atoms):
			line = np.mean(self.slopes[index])*np.asarray(self.timesteps)+np.mean(self.intercepts[index])
			plt.plot(self.timesteps, line, color='C%d'%(index), label='Mean : %s'%(self.types_of_atoms[index]))
			if self.data == True:
				print('---')
				print('Mean Diffusion Coefficient (X, Y and Z) : %s = %.10f cm^2/s, %.10f m^2/s; Standard Deviation = %.10f cm^2/s, %.10f m^2/s' % (self.types_of_atoms[index],np.mean(self.slopes[index])*(0.1), np.mean(self.slopes[index])*(10**-5), np.std(self.slopes[index])*(0.1), np.std(self.slopes[index])*(10**-5)))
				print('---')




		# This plots the lines between the segments
		for seg in self.segment_line:
			plt.plot([seg,seg],[np.amin(self.cont_xyz_segment_ensemble_average), np.amax(self.cont_xyz_segment_ensemble_average)], color='grey')
		if self.molecule_index == None:
			line = np.mean(self.slopes)*np.asarray(self.timesteps)+np.mean(self.intercepts)
			plt.plot(self.timesteps, line, color='black',label='Mean')
			diff=np.mean(self.slopes)*(10**-5)
			if self.data == True:
				print('---')
				print('Mean Diffusion Coefficient (X, Y and Z) = %.10f cm^2/s, %.10f m^2/s; Standard Deviation = %.10f cm^2/s, %.10f m^2/s' % (np.mean(self.slopes)*(0.1), diff, np.std(self.slopes)*(0.1), np.std(self.slopes)*(10**-5)))
				print('---')

		plt.legend(loc='best')
		plt.xlabel('Time (fs)')
		plt.ylabel(r'Mean Square Displacement ($\AA^2$)')

		if self.plot_graph == True:
			plt.show()
		return diff
