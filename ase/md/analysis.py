from math import floor
import numpy as np

class DiffusionCoefficient:

	def __init__(self, traj, timestep, steps_between_saved_images, atom_index=None, data=False):
		'''
		Calculates Diffusion Coefficients for atoms and molecules

		Parameters:
			traj (Trajectory): Trajectory of atoms objects (images)
			timestep (Int): Timestep used between each images
			steps_between_saved_images (Int): Interval used when writing the .traj file 
			molecule_index (List of Int): The indices of atoms whose Diffusion Coefficient is to be calculated
			data (Bool): If the user wants to print the data or not

			This function calculates the Diffusion Coefficient for the given .traj file using the Einstein Equation:
			⟨|r(t)−r(0)|**2⟩=6Dt (where r(t) is the position of atom at time t, D is the Diffusion Coefficient)
			wiki : https://en.wikibooks.org/wiki/Molecular_Simulation/Diffusion_Coefficients
		'''

		self.traj = traj
		self.timestep = timestep
		self.steps_between_saved_images = steps_between_saved_images
		self.atom_index = atom_index
		self.data = data

		# Condition used if user wants to calculate Diffusion Coefficients for a specific atom or group of atoms (molecule)
		if self.atom_index == None:
			self.types_of_atoms = list(set(traj[0].get_chemical_symbols()))
			self.no_of_types_of_atoms = len(self.types_of_atoms)
			self.no_of_atoms = [traj[0].get_chemical_symbols().count(symbol) for symbol in self.types_of_atoms]
		else:
			# This isn't strictly ideal - it only deals with one atom in a system. Better to generalise this long-term
			self.no_of_types_of_atoms = 1
			self.no_of_atoms = 1
			self.types_of_atoms = [traj[0].symbols[self.atom_index]]

	def initialise_lists(self):

		self.time_between_images = self.timestep * self.steps_between_saved_images
		self.total_images = len(self.traj) - self.ignore_n_images
		# List of Timestep used for plotting (x-axis)
		self.timesteps = list(np.arange(0,self.total_images*self.time_between_images,self.time_between_images))

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

		self.slopes=[]
		self.intercepts=[]
		for i in range(self.no_of_types_of_atoms):
			self.slopes.append([])
			self.intercepts.append([])

		self.cont_xyz_segment_ensemble_average = 0
		self.segment_line=[]

	def calculate(self, ignore_n_images=0, number_of_segments=1):
		'''
		Calculate the diffusion coefficients.

		Parameter:
			ignore_n_images (Int): Number of images you want to ignore from the start of the trajectory, e.g. during equilibration
			no_of_segments (Int): Divides the given trajectory in to segments to allow statistical analysis
		'''

		self.ignore_n_images = ignore_n_images
		self.no_of_segments = number_of_segments

		self.initialise_lists()

		for segment_no in range(self.no_of_segments):
			start = self.ignore_n_images + (segment_no*floor(self.total_images/self.no_of_segments))  
			end = start + floor(self.total_images/self.no_of_segments)
			seg = self.traj[start:end]

			# For each image
			for image_no in range(1,len(seg)): 
				xyz_disp = [[0., 0., 0.]*self.no_of_types_of_atoms]
				
				# Calculating for all the atoms
				if self.atom_index == None: 
					# For each atom, work out displacement from start coordinate and collect information with like atoms
					for atom_no in range(len(seg[image_no])):
						sym_index = self.types_of_atoms.index(seg[image_no].symbols[atom_no])
						for xyz in range(3):
							xyz_disp[sym_index][xyz] += np.square(seg[image_no].positions[atom_no][xyz] - seg[0].positions[atom_no][xyz])
					
					# For each atom species, use xyz_disp to calculate the average data
					for sym_index in range(self.no_of_types_of_atoms):
						# This is the average displacement for an atom species over the entire segment
						self.segment_ensemble_average[segment_no][sym_index].append(np.sum(xyz_disp[sym_index])/(6*self.no_of_atoms[sym_index]))
						# This is the average displacement in X, Y and Z over the entire segment
						for xyz in range(3):
							self.xyz_segment_ensemble_average[segment_no][sym_index][xyz].append(xyz_disp[sym_index][xyz]/(2*self.no_of_atoms[sym_index]))

						

		
				else: # Calculating for specific atom or group of atoms (molecule)
					com_orig = np.zeros(3)
					com_new = np.zeros(3)
					for atom_no in self.atom_index:
						for xyz in range(3):
							com_orig[xyz] += seg[0].positions[atom_no][xyz]
							com_new[xyz] += seg[image_no].positions[atom_no][xyz]
					com_orig /= len(self.atom_index)
					com_new /= len(self.atom_index)
					for xyz in range(3):
						xyz_disp[0][xyz] += np.square(com_new[xyz] - com_orig[xyz])
						self.xyz_segment_ensemble_average[segment_no][0][xyz].append(xyz_disp[0][xyz]/2)
					self.segment_ensemble_average[segment_no][0].append(np.sum(xyz_disp[0][xyz]/6))

		self.cont_xyz_segment_ensemble_average = self.xyz_segment_ensemble_average

		for segment_no in range(self.no_of_segments):
			if segment_no>0:
			    for index in range(self.no_of_types_of_atoms):
			        for xyz in range(3):
			            self.cont_xyz_segment_ensemble_average[segment_no][index][xyz]=np.array(self.cont_xyz_segment_ensemble_average[segment_no][index][xyz])+self.cont_xyz_segment_ensemble_average[segment_no-1][index][xyz][-1]
			    
			for index in range(self.no_of_types_of_atoms):	
				for xyz in range(3):
					ratio = floor(self.total_images/self.no_of_segments)
					slope, intercept = self.fit_data(self.timesteps[segment_no*ratio+1:segment_no*ratio+ratio], self.xyz_segment_ensemble_average[segment_no][index][xyz])
					self.slopes[index].append(slope)
					self.intercepts[index].append(intercept)

			self.segment_line.append(self.timesteps[segment_no*floor(self.total_images/self.no_of_segments)+floor(self.total_images/self.no_of_segments)-1] + self.time_between_images)

		return np.mean(self.slopes)*(10**-5)

	def fit_data(self,x, y):
		# Moved local to usage
		from scipy.stats import linregress

		# Generated linear fit  
		slope, intercept, r_value, p_value, std_err = linregress(x,y)

		return slope, intercept

	def plot(self, continuous=True):
		'''
		Auto-plot of Diffusion Coefficient data

		Parameters:
			continuous (Bool): Set True if the user wants the graph to be continuous
		'''

		# Moved matplotlib into the function so it is not loaded unless needed
		import matplotlib.pyplot as plt
		
		color_list = plt.cm.Set3(np.linspace(0, 1, self.no_of_types_of_atoms))
		xyz_labels=['X','Y','Z']

		# We need to separate things out here so all the calculations are done in "calculate", and then the stored data is plotted with "plot"
		# Calculate should return the diffusion coefficient, not plot.
		if len(self.slopes) == 0:
			diff = self.calculate()
		
		for segment_no in range(self.no_of_segments):
		    
			for index in range(self.no_of_types_of_atoms):	
				for xyz in range(3):
					if segment_no == 0:
						custom_label = 'Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], xyz_labels[xyz])
					else:
						custom_label = '_Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], xyz_labels[xyz]) # To remove label for other segments as all the segments are of same colour
					ratio = floor(self.total_images/self.no_of_segments)
					self.plot_data(self.timesteps[segment_no*ratio+1:segment_no*ratio+ratio], self.xyz_segment_ensemble_average[segment_no][index][xyz], self.slopes[index], self.intercept[index], plt, color_list[index], custom_label, continuous, self.cont_xyz_segment_ensemble_average[segment_no][index][xyz])

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
		if self.atom_index == None:
			line = np.mean(self.slopes)*np.asarray(self.timesteps)+np.mean(self.intercepts)
			diff=np.mean(self.slopes)*(10**-5)
			plt.plot(self.timesteps, line, color='black',label='Mean')
			if self.data == True:
				print('---')
				print('Mean Diffusion Coefficient (X, Y and Z) = %.10f cm^2/s, %.10f m^2/s; Standard Deviation = %.10f cm^2/s, %.10f m^2/s' % (np.mean(self.slopes)*(0.1), diff, np.std(self.slopes)*(0.1), np.std(self.slopes)*(10**-5)))
				print('---')

		plt.legend(loc='best')
		plt.xlabel('Time (fs)')
		plt.ylabel(r'Mean Square Displacement ($\AA^2$)')

		plt.show()

	def plot_data(self, x, y, slope, intercept, plt, color, label=None, cont=False, cont_y=None):
		# Moved this local to its use.
		xyz_markers = {'X':'o', 'Y':'s','Z':'^'}

		#line = slope*np.asarray(x)+intercept
		if label[0] != '_':
			skip = label.index('A')
		else:
			skip = 0

		mark = xyz_markers[label[-1]]

		if cont == False:
			plt.plot(x, y, color=color, marker=mark,label=label[skip:], linewidth=0)
		else:
			plt.plot(x, cont_y, color=color, marker=mark,label=label[skip:], linewidth=0)

		if label[0] == '_':
			label=label[1:]

		# print(r'Intercept = %.10f $\AA^2$; Diffusion Coefficient = %.10f $\AA^2$/fs' % (intercept, slope))
		# Converting to more common units of cm^2/s => multiply by 10^-1
		# \AA^2 => cm^2 requires multiplying by (10^-8)^-2
		# fs => s requires dividing by 10^-15
		if self.data == True:
			print('---')
			print(r'%10s: Intercept = %.10f cm^2; Diffusion Coefficient = %.10f cm^2/s, %.10f m^2/s' % (label, intercept/(10**8), slope*(0.1), slope*(10**-5)))
