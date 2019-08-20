import numpy as np

class DiffusionCoefficient:

    def __init__(self, traj, timestep, steps_between_saved_images, atom_indices=None, molecule=False):
        '''
        Calculates Diffusion Coefficients for atoms and molecules

        Parameters:
            traj (Trajectory): 
                Trajectory of atoms objects (images)
            timestep (Int): 
                Timestep used between each images
            steps_between_saved_images (Int): 
                Interval used when writing the .traj file 
            molecule_index (List of Int): 
                The indices of atoms whose Diffusion Coefficient is to be calculated

            This function calculates the Diffusion Coefficient for the given .traj file using the Einstein Equation:
            ⟨|r(t)−r(0)|**2⟩ = 6Dt (where r(t) is the position of atom at time t, D is the Diffusion Coefficient)
            Solved herein as y = mx + c, i.e. 1/6 ⟨|r(t)−r(0)|**2⟩ = Dt, so m = D, c = 0

            wiki : https://en.wikibooks.org/wiki/Molecular_Simulation/Diffusion_Coefficients
        '''

        self.traj = traj
        self.timestep = timestep
        self.steps_between_saved_images = steps_between_saved_images

        # Condition used if user wants to calculate diffusion coefficients for specific atoms or all atoms
        self.atom_indices = atom_indices
        if self.atom_indices == None:
            self.atom_indices = [i for i in range(len(traj[0]))] 

        # Condition if we are working with the mobility of a molecule, need to manage arrays slightly differently
        self.molecule = molecule
        if molecule:
            self.types_of_atoms = ["molecule"]
            self.no_of_atoms = [1]
        else:
            self.types_of_atoms = list(set(traj[0].get_chemical_symbols()))
            self.no_of_atoms = [traj[0].get_chemical_symbols().count(symbol) for symbol in self.types_of_atoms]

        self.no_of_types_of_atoms = len(self.types_of_atoms)

    def initialise_arrays(self, ignore_n_images, number_of_segments):

        from math import floor
        total_images = len(self.traj) - ignore_n_images
        self.no_of_segments = number_of_segments
        self.len_segments = floor(total_images/self.no_of_segments)

        time_between_images = self.timestep * self.steps_between_saved_images
        # These are the data objects we need when plotting information. First the x-axis, timesteps
        self.timesteps = np.arange(0,total_images*time_between_images,time_between_images)
        # This holds all the data points for the diffusion coefficients, averaged over atoms
        self.xyz_segment_ensemble_average = np.zeros((self.no_of_segments,self.no_of_types_of_atoms,3,self.len_segments))
        # This holds all the information on linear fits, from which we get the diffusion coefficients
        self.slopes = np.zeros((self.no_of_segments,self.no_of_types_of_atoms,3))
        self.intercepts = np.zeros((self.no_of_segments,self.no_of_types_of_atoms,3))

        self.cont_xyz_segment_ensemble_average = 0

    def calculate(self, ignore_n_images=0, number_of_segments=1):
        '''
        Calculate the diffusion coefficients.

        Parameter:
            ignore_n_images (Int): 
                Number of images you want to ignore from the start of the trajectory, e.g. during equilibration
            number_of_segments (Int): 
                Divides the given trajectory in to segments to allow statistical analysis
        '''

        # Setup all the arrays we need to store information
        self.initialise_arrays(ignore_n_images, number_of_segments)

        for segment_no in range(self.no_of_segments):
            start = segment_no*self.len_segments  
            end = start + self.len_segments
            seg = self.traj[ignore_n_images+start:ignore_n_images+end]

            # If we are considering a molecular system, work out the COM for the starting structure
            if self.molecule:
                com_orig = np.zeros(3)
                for atom_no in self.atom_indices:
                    com_orig[:] += seg[0].positions[atom_no][:] / len(self.atom_indices)

            # For each image, calculate displacement.
            # I spent some time deciding if this should run from 0 or 1, as the displacement will be zero for 
            # t = 0, but this is a data point that needs fitting too and so should be included
            for image_no in range(0,len(seg)): 
                # This object collects the xyz displacements for all atom species in the image
                xyz_disp = [[0., 0., 0.]*self.no_of_types_of_atoms]
                
                # Calculating for each atom individually, grouping by species type (e.g. solid state)
                if not self.molecule:
                    # For each atom, work out displacement from start coordinate and collect information with like atoms
                    for atom_no in self.atom_indices:
                        sym_index = self.types_of_atoms.index(seg[image_no].symbols[atom_no])
                        xyz_disp[sym_index][:] += np.square(seg[image_no].positions[atom_no][:] - seg[0].positions[atom_no][:])
        
                else: # Calculating for group of atoms (molecule) and work out squared displacement
                    com_disp = np.zeros(3)
                    for atom_no in self.atom_indices:
                        com_disp[:] += seg[image_no].positions[atom_no][:] / len(self.atom_indices)
                    xyz_disp[0][:] += np.square(com_disp[:] - com_orig[:])

                # For each atom species or molecule, use xyz_disp to calculate the average data                      
                for sym_index in range(self.no_of_types_of_atoms):
                    # This is the average displacement in X, Y and Z for each species in the entire segment
 		    # Normalise by degrees of freedom and number of atoms.                            
                    denominator = (2*self.no_of_atoms[sym_index])
                    for xyz in range(3):
                        self.xyz_segment_ensemble_average[segment_no][sym_index][xyz][image_no] = (xyz_disp[sym_index][xyz]/denominator)

            # We've collected all the data for this entire segment, so now to fit the data.
            for sym_index in range(self.no_of_types_of_atoms):    
                self.slopes[segment_no][sym_index], self.intercepts[segment_no][sym_index] = self.fit_data(self.timesteps[start:end], self.xyz_segment_ensemble_average[segment_no][sym_index][:])

        self.cont_xyz_segment_ensemble_average = self.xyz_segment_ensemble_average.copy()

        #return self.slopes
        return np.mean(self.slopes)*(10**-5)

    def fit_data(self, x, y):

        # Simpler implementation but disabled as fails Conda tests.
        # from scipy.stats import linregress
        # slope, intercept, r_value, p_value, std_err = linregress(x,y)
       
        # Initialise objects
        slopes = np.zeros(3)
        intercepts = np.zeros(3)

        # Convert into suitable format for lstsq
        x_edited = np.vstack([np.array(x), np.ones(len(x))]).T
        # Calculate slopes for x, y and z-axes
        for xyz in range(3):
            slopes[xyz], intercepts[xyz] = np.linalg.lstsq(x_edited, y[xyz], rcond=-1)[0]

        return slopes, intercepts

    def plot(self, continuous=True, print_data=False):
        '''
        Auto-plot of Diffusion Coefficient data

        Parameters:
            continuous (Bool): 
                Set True if the user wants the graph to be continuous
            print_data (Bool): 
                Set True to get output information for details of diffusion coefficient plots
        '''

        # Moved matplotlib into the function so it is not loaded unless needed
        # Could be provided as an input variable, so user can work with it further?
        import matplotlib.pyplot as plt
        
        # Define some aesthetic variables
        color_list = plt.cm.Set3(np.linspace(0, 1, self.no_of_types_of_atoms))
        xyz_labels=['X','Y','Z']
        xyz_markers = {'X':'o', 'Y':'s','Z':'^'}

        # Check if we have data to plot, if not calculate it.
        if len(self.slopes) == 0:
            diff = self.calculate()
        
	# AL: Still tidying but we are nearly there.
        for segment_no in range(self.no_of_segments):
            start = segment_no*self.len_segments  
            end = start + self.len_segments
            
            for index in range(self.no_of_types_of_atoms):    
                for xyz in range(3):
                    if segment_no == 0:
                        custom_label = 'Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], xyz_labels[xyz])
                    else:
                        custom_label = '_Segment No. - %d ; Atom - %s ; %s'%(segment_no+1, self.types_of_atoms[index], xyz_labels[xyz]) # To remove label for other segments as all the segments are of same colour

                    if custom_label[0] != '_':
                        skip = custom_label.index('A')
                    else:
                        skip = 0

                    mark = xyz_markers[custom_label[-1]]

                    if continuous == False:
                       plt.plot(self.timesteps[start:end], self.xyz_segment_ensemble_average[segment_no][index][xyz], color=color_list[index], marker=mark, label=custom_label[skip:], linewidth=0)
                    else:
                       plt.plot(self.timesteps[start:end], self.cont_xyz_segment_ensemble_average[segment_no][index][xyz], color=color_list[index], marker=mark, label=custom_label[skip:], linewidth=0)

                    if custom_label[0] == '_':
                        custom_label=custom_label[1:]

                    # print(r'Intercept = %.10f $\AA^2$; Diffusion Coefficient = %.10f $\AA^2$/fs' % (intercept, slope))
                    # Converting to more common units of cm^2/s => multiply by 10^-1
                    # \AA^2 => cm^2 requires multiplying by (10^-8)^-2
                    # fs => s requires dividing by 10^-15
                    if print_data:
                        print('---')
                        print(r'%10s: Intercept = %.10f cm^2; Diffusion Coefficient = %.10f cm^2/s, %.10f m^2/s' % (label, intercept/(10**8), slope*(0.1), slope*(10**-5)))

        for index in range(self.no_of_types_of_atoms):
            line = np.mean(self.slopes[index])*self.timesteps+np.mean(self.intercepts[index])
            plt.plot(self.timesteps, line, color='C%d'%(index), label='Mean : %s'%(self.types_of_atoms[index]))
            if print_data:
                print('---')
                print('Mean Diffusion Coefficient (X, Y and Z) : %s = %.10f cm^2/s, %.10f m^2/s; Standard Deviation = %.10f cm^2/s, %.10f m^2/s' % (self.types_of_atoms[index],np.mean(self.slopes[index])*(0.1), np.mean(self.slopes[index])*(10**-5), np.std(self.slopes[index])*(0.1), np.std(self.slopes[index])*(10**-5)))
                print('---')

        # This plots the lines between the segments
        for segment_no in range(self.no_of_segments):
            x_coordinate = self.timesteps[segment_no*self.len_segments+self.len_segments-1] + self.timestep * self.steps_between_saved_images
            plt.plot([x_coordinate, x_coordinate],[np.amin(self.cont_xyz_segment_ensemble_average), np.amax(self.cont_xyz_segment_ensemble_average)], color='grey')
        if self.atom_indices == None:
            line = np.mean(self.slopes)*self.timesteps+np.mean(self.intercepts)
            diff=np.mean(self.slopes)*(10**-5)
            plt.plot(self.timesteps, line, color='black',label='Mean')
            if print_data:
                print('---')
                print('Mean Diffusion Coefficient (X, Y and Z) = %.10f cm^2/s, %.10f m^2/s; Standard Deviation = %.10f cm^2/s, %.10f m^2/s' % (np.mean(self.slopes)*(0.1), diff, np.std(self.slopes)*(0.1), np.std(self.slopes)*(10**-5)))
                print('---')

        plt.legend(loc='best')
        plt.xlabel('Time (fs)')
        plt.ylabel(r'Mean Square Displacement ($\AA^2$)')

        plt.show()

    def plot_data(self, x, y, plt, color, label=None, cont=False, cont_y=None):

        if custom_label[0] != '_':
            skip = custom_label.index('A')
        else:
            skip = 0

        mark = xyz_markers[custom_label[-1]]

        if cont == False:
            plt.plot(x, y, color=color, marker=mark,label=label[skip:], linewidth=0)
        else:
            plt.plot(x, cont_y, color=color, marker=mark,label=label[skip:], linewidth=0)

        if label[0] == '_':
            label=label[1:]

        return label