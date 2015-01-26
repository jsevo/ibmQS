package nl.tudelft.bt.model;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

import nl.tudelft.bt.model.apps.output.SimulationResultsWriter;
import nl.tudelft.bt.model.apps.output.SolidsConcentrationWriter;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.detachment.levelset.DetachedBiomassContainer;
import nl.tudelft.bt.model.exceptions.BiofilmProblemException;
import nl.tudelft.bt.model.exceptions.InvalidValueException;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;
import nl.tudelft.bt.model.multigrid.BiomassBooleanMatrix;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;
import nl.tudelft.bt.model.parsers.SimulationResultParser;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.reaction.Reaction;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;
import nl.tudelft.bt.model.util.UnlimitedFloatArray;

/**
 * Main model class. Implements singleton design pattern to impose a unique
 * instance of model. 22042004 new change
 * 
 * @author Jo?o Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Model implements Serializable {
	// the singleton instance of model (exists 1 and only 1)
	private static Model _model = new Model();

	private int _seed; // random number generator seed

	private Random _rng; // random number generator

	// computation
	private int _dimensionality;

	private int _iterationCounter;

	private float _time;

	private long _realInitialTime;

	private float _timeStep;

	private float _compulsoryTimeStep;

	private float _maximumTimeStep;

	private float _minimumTimeStep;

	private boolean _writeTimedWriters; // flag for timed writers

	private boolean _overrideTimeStep;
	
	//added by Jonas
	private boolean _overrideTimeStepWithSpecificTimesForAllData;
	private int _tOut;
	//end added by Jonas
	private ArrayList<TimeStepConstraint> _timeStepConstraints;

	// attributes for simulation finishing criterium
	private float _finishIterationTime;

	private float _finalBiofilmHeight;

	private float _maxBiovolume;

	private float _maxRunLength;

	private boolean _finishSimulation;

	private UnlimitedFloatArray _timeSeries;

	public ContinuousCoordinate systemSize;

	public float referenceSystemSide;

	private float _systemGridSide;

	private int _numberOfGridNodes;

	public BiomassContainer biomassContainer;

	private BiomassSpecies[] _biomassSpecies;

	private SoluteSpecies[] _soluteSpecies;
	
	// Jonas
	private SoluteSpecies _nutrientForConsumptionStopCriterion;
	private float _maxConsumption;
	private float _cummulativeConsumptionOfNutrient;
	// end Jonas

	private ParticulateSpecies[] _particulateSpecies;

	private ArrayList<Reaction> _reactions;

	private BoundaryLayer _boundaryLayer;

	private float _maxBiofilmHeight;

	// detached biomass
	private DetachedBiomassContainer _detachedBiomass;

	// global system scale
	private float _carrierArea;

	private float _residenceTime;

	private float _reactorVolume;

	/**
	 * Model constructor (private to implement Singleton Model), initializes
	 * random seed to its default value.
	 */
	private Model() {
		// initialize the random number generator
		// to default seed value
		// predifined
		// _seed = 1007;
		// random
		_seed = (int) System.currentTimeMillis();
		reset();
	}

	/**
	 * Returns the Singleton instance of model
	 * 
	 * @return the singleton instace of model
	 */
	public static Model model() {
		return _model;
	}

	/**
	 * Resets the model.
	 */
	public void reset() {
		// initialize the random number generator
		_rng = new java.util.Random(_seed); // use this for fixed seed
		// TODO delete:
		System.out.println("seed " + _seed);
		// _rng = new java.util.Random(); // use this for variable seed
		_dimensionality = 0;
		_iterationCounter = 0;
		_realInitialTime = System.currentTimeMillis(); // save the initial time
		// for benchmark purposes
		_time = 0;
		_timeStep = 0;
		_compulsoryTimeStep = 0;
		_maximumTimeStep = 0;
		_minimumTimeStep = 0;
		_writeTimedWriters = true;
		_overrideTimeStep = false;
		// added by Jonas
		_overrideTimeStepWithSpecificTimesForAllData = false;
		_tOut = 0;
		_maxConsumption = Float.POSITIVE_INFINITY;
		_nutrientForConsumptionStopCriterion = null;
		_cummulativeConsumptionOfNutrient = 0f;
		// end added by Jonas
		_finishIterationTime = Float.POSITIVE_INFINITY;
		_maxBiovolume = Float.POSITIVE_INFINITY;
		_maxRunLength = Float.POSITIVE_INFINITY;
		_finalBiofilmHeight = Float.POSITIVE_INFINITY;
		_finishSimulation = false;
		_timeSeries = new UnlimitedFloatArray();
		_timeSeries.add(_time);
		referenceSystemSide = 0;
		_maxBiofilmHeight = 0;
		systemSize = null;
		_timeStepConstraints = new ArrayList<TimeStepConstraint>();
		// delete all references to bacteria and the shoving grid:
		biomassContainer = null;
		// removes all discrete multigrid data on species
		_biomassSpecies = null;
		_soluteSpecies = null;
		// the biofilm reactor propetries
		_residenceTime = 0;
		_carrierArea = 0;
		_reactorVolume = 0;
		_detachedBiomass = new DetachedBiomassContainer(_timeSeries);
		// Initialize the reactions container
		_reactions = new ArrayList<Reaction>();
	}

	/**
	 * Set the value for the random number generator seed
	 * 
	 * @param i
	 *            value to set seed to
	 */
	public void setSeed(int i) {
		_seed = i;
	}

	/**
	 * Get the value for the random number generator seed
	 * 
	 * @param i
	 *            value to set seed to
	 */
	public int getSeed() {
		return _seed;
	}

	/**
	 * The only method o provides a random number that should be used anywhere
	 * where random number generation is require in the application.
	 * 
	 * @return random float between 0 and 1
	 * @see Random
	 */
	public float getRandom() {
		return _rng.nextFloat();
	}

	/**
	 * The only method o provides a random number that should be used to
	 * generate a normally distributed random variavle.
	 * 
	 * @return random with average 0 and standard deviation 1
	 * @see Random
	 */
	public float getRandomFromNormalDistribution() {
		return (float) (_rng.nextGaussian());
	}

	/**
	 * Builds the spatial system
	 * 
	 * @param d
	 *            dimensionality (2D or 3D)
	 * @param s
	 *            system size
	 * @param g
	 *            grid size for solute concentration field
	 * @throws InvalidValueException
	 *             id ((d != 2) && (d != 3)) and if g is not of form 2^n + 1
	 */
	public void buildSystem(int d, float s, int g) throws InvalidValueException {
		if ((d != 2) && (d != 3))
			throw new InvalidValueException("dimensionality value (" + d
					+ ") not valid");
		_dimensionality = d;
		systemSize = new ContinuousCoordinate(s, s, (d == 3 ? s : 0));
		referenceSystemSide = s;
		_systemGridSide = s / g;
		_numberOfGridNodes = g;
		_maxBiofilmHeight = s;
		// define the grid resolution
		MultigridVariable.setGrid(g, g, (d == 3 ? g : 1));
	}

	/**
	 * Get value of time step
	 * 
	 * @return time step
	 */
	public float getTimeStep() {
		return _timeStep;
	}

	/**
	 * Set the size for biofilm detachment (thickness if flate geometry biofilm,
	 * radius if granule or colony
	 * 
	 * @param h
	 *            maximum biofilm size (units of length)
	 */
	public void setMaximumBiofilmSize(float h) {
		if (h <= _model.systemSize.x)
			_maxBiofilmHeight = h;
		else
			throw (new InvalidValueException(
					"Proposed value for maximum  height (" + h
							+ ") is higher than system height ("
							+ _model.systemSize.x + ")"));
	}

	/**
	 * Get the biofilm detachment height
	 * 
	 * @return height biofilm detachment height
	 */
	public float getMaximumBiofilmHeight() {
		return _maxBiofilmHeight;
	}

	/**
	 * Get the current biofilm height
	 * 
	 * @return current biofilm height
	 */
	public float getCurrentBiofilmHeight() {
		return biomassContainer.getTopBiofilmHeight();
	}

	/**
	 * Get the current biovolume. Biovolume is computed from the multigrid
	 * 
	 * @return current biovolume
	 */
	public float getCurrentBiovolume() {
		return biomassContainer.computeBiovolume();
		// return MultigridVariable.computeBiovolume(_particulateSpecies);
	}

	/**
	 * @return the total fixed biomass [kg/m2]
	 */
	public float getCurrentTotalBiomass() {
		return biomassContainer.getTotalBiomass();
	}

	/**
	 * Set the array of bacteria species, chemical species and the boundary
	 * layer object that will be used by the multigrid solver.
	 * 
	 * @param b
	 *            array of bacteria species
	 * @param c
	 *            array of bacteria species
	 * @param bl
	 *            boundary layer object
	 */
	public void setupDiffusionReactionSystem(BiomassSpecies[] b,
			SoluteSpecies[] c, BoundaryLayer bl) {
		_biomassSpecies = b;
		_soluteSpecies = c;
		_boundaryLayer = bl;
		// get the array of fixed substances from the biomass species array
		Collection<ParticulateSpecies> col = new ArrayList<ParticulateSpecies>();
		for (int i = 0; i < _biomassSpecies.length; i++) {
			Collection<ParticulateSpecies> clocal = _biomassSpecies[i]
					.getFixedSpeciesAsCollection();
			col.removeAll(clocal);
			col.addAll(clocal);
		}
		_particulateSpecies = new ParticulateSpecies[col.size()];
		int i = 0;
		for (Iterator<ParticulateSpecies> iter = col.iterator(); iter.hasNext();) {
			_particulateSpecies[i++] = (ParticulateSpecies) iter.next();
		}
	}

	/**
	 * Add a reaction to the reaction container. This method is called directly
	 * by the Reaction contructor, so that each reaction created is
	 * automatically stored
	 * 
	 * @param r
	 *            reactioin to add to model
	 */
	public void addReaction(Reaction r) {
		_reactions.add(r);
	}

	/**
	 * Creates the bacteria container as a spherical particle container.
	 * 
	 * @param maximumRadius
	 *            maximum radius of biomass particles (for division)
	 * @param minimumRadius
	 *            minimum radius of biomass particles (for removal)
	 * @param kshov
	 *            shoving multiplier
	 * @param fshov
	 *            fraction of remaining particles
	 * @param rdetach
	 *            detachment rate [10^-15g/um^2/h]
	 */
	public void buildBacteriaContainer(float maximumRadius,
			float minimumRadius, float kshov, float fshov) {
		biomassContainer = new BiomassParticleContainer(maximumRadius,
				minimumRadius, kshov, fshov);
	}




	/**
	 * Inoculate system at random positions with arbitrary number of cells of
	 * each defined species. Mass of cell is the average mass defined in species
	 * 
	 * @param n
	 *            number of cells to inoculate for each species
	 * @throws NonMatchingNumberException
	 *             if length of n does not match number of defined bacteria
	 *             species
	 */
	public void inoculateRandomlyMultispecies(int[] n)
			throws NonMatchingNumberException {
		if (n.length != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// Inoculate cells randomly (shhuffling order)
        int nspecies = n.length;
        int[] numberOfCellsToPlace = new int[n.length];
        int totalNumberOfCellsToPlace = 0;
        float[] ratioOfCellsToPlace = new float[n.length];
        for (int i = 0; i < nspecies; i++) {
                numberOfCellsToPlace[i] = n[i];
                totalNumberOfCellsToPlace += n[i];
        }	
		while (true) {
			// get the random value
			// PREVIOUS CODE: BUG! Should depend on frequency of
            // each species
            //int r = (int) Math.floor(getRandom() * nspecies);
			//
			//Jonas' addition to Sara's bug fix below:
			//Recalculation of the probability to seed a certain species based on what has been seeded thus far.
			totalNumberOfCellsToPlace = 0;
	        for (int i = 0; i < nspecies; i++) {
                totalNumberOfCellsToPlace += numberOfCellsToPlace[i];} //numberOfCellsToPlace gets reduced whenever a cell is seeded.
			for (int i = 0; i < nspecies; i++){
                ratioOfCellsToPlace[i] = numberOfCellsToPlace[i]/(float)totalNumberOfCellsToPlace;} //new probabilities
			//End of Jonas' addition to Sara's bug fix below
            //
			// NEW IMPLEMENTATION (by Sara):
            // Frequency-dependent inoculation order
            float rand = getRandom();
            float ratioCounter = 0;
            int r = -1;
            for (int i = 0; i < nspecies; i++) {
                    ratioCounter += ratioOfCellsToPlace[i];
                    if (rand < ratioCounter) { r = i;break; }
            }
            // end of NEW IMPLEMENTATION (by Sara)
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float y = getRandom() * systemSize.y;
					float z = getRandom() * systemSize.z;
					// add to the bacteria list
					biomassContainer
							.placeBiomassAt(_biomassSpecies[r], 0, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}
		spreadByShovingCompletely();
	}
	
	
	
	
	
	/**
	 * Inoculate system at random positions with arbitrary number of cells of
	 * each defined species. Mass of cell is the average mass defined in species
	 * @author Jonas
	 * @param n
	 *            number of cells to inoculate for each species
	 * @param h minimal x location (height above surface)
	 * @param y1,y2 these are the left and right borders of the region where to inoculate
	 * @throws NonMatchingNumberException
	 *             if length of n does not match number of defined bacteria
	 *             species
	 */
	public void inoculateRandomlyMultispeciesAtMinHeightBetweenY1andY2(int[] n, float h,float y1, float y2)
			throws NonMatchingNumberException {
		if (n.length != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// Inoculate cells randomly (shhuffling order)
		int nspecies = n.length;
		int[] numberOfCellsToPlace = new int[n.length];
		int totalNumberOfCellsToPlace = 0;
		float[] ratioOfCellsToPlace = new float[n.length];
		for (int i = 0; i < nspecies; i++) {
                numberOfCellsToPlace[i] = n[i];
                totalNumberOfCellsToPlace += n[i];
        }	
		while (true) {
			// get the random value
			// PREVIOUS CODE: BUG! Should depend on frequency of
            // each species
            //int r = (int) Math.floor(getRandom() * nspecies);
			//
			//Jonas' addition to Sara's bug fix below:
			//Recalculation of the probability to seed a certain species based on what has been seeded thus far.
			totalNumberOfCellsToPlace = 0;
	        for (int i = 0; i < nspecies; i++) {
                totalNumberOfCellsToPlace += numberOfCellsToPlace[i];} //numberOfCellsToPlace gets reduced whenever a cell is seeded.
			for (int i = 0; i < nspecies; i++){
                ratioOfCellsToPlace[i] = numberOfCellsToPlace[i]/(float)totalNumberOfCellsToPlace;} //new probabilities
			//End of Jonas' addition to Sara's bug fix below
            //
			// NEW IMPLEMENTATION (by Sara):
            // Frequency-dependent inoculation order
            float rand = getRandom();
            float ratioCounter = 0;
            int r = -1;
            for (int i = 0; i < nspecies; i++) {
                    ratioCounter += ratioOfCellsToPlace[i];
                    if (rand < ratioCounter) { r = i;break; }
            }
            // end of NEW IMPLEMENTATION (by Sara)
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float y = getRandom() * (y2-y1)+y1;
					float z = getRandom() * systemSize.z;
					// add to the bacteria list
					biomassContainer
							.placeBiomassAt(_biomassSpecies[r], h, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done){
				System.out.println(numberOfCellsToPlace[0]);
				break;
			}
		}
		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	

	/**
	 * Place a single cell in center of system with masses defined by the
	 * default array of masses of the biomass species
	 * 
	 * @param s
	 *            species of biomass to add
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @return reference to particle added
	 */
	public BiomassParticle placeSingleCellInCenter(BiomassSpecies s) {
		float x = 0;
		float y = systemSize.y * 0.5f;
		float z = systemSize.z * 0.5f;
		return biomassContainer.placeBiomassAt(s, x, y, z);
	}

	/**
	 * Place biomass in the system with masses defined by the default array of
	 * masses of the biomass species
	 * 
	 * @param s
	 *            species of biomass to add
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @return reference to particle added
	 */
	public BiomassParticle placeBiomass(BiomassSpecies s, float x, float y,
			float z) {
		return biomassContainer.placeBiomassAt(s, x, y, z);
	}





	/**
	 * Place randonmly cells of each species everythere in space
	 * 
	 * @param n
	 * @param rad
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesEverywhere(int[] n, float rad)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float x = getRandom() * systemSize.x;
					float y = getRandom() * systemSize.y;
					float z = getRandom() * systemSize.z;
					// add to the bacteria list
					biomassContainer
							.placeBiomassAt(_biomassSpecies[r], x, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();

	}


	/**
	 * Update the concentration of biomass in the diffusion-reaction data
	 * matrices.
	 */
	public void updateBioDiscreteData() {
		// reset values of all BacteriaSpecies
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].resetDiscreteMatrix();
		}
		// iterate through array of bacteria and increase
		// concentration value in corresponding BacteriaSpecies
		biomassContainer.addContributionsToSpeciesDiscreteData();
		// update boundary values of all BacteriaSpecies
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].refreshBoundaryConditions();
		}
	}

	/**
	 * Compute the initial concentration fields and initializaes the global
	 * rates time series
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void initializeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		updateBioDiscreteData();
		MultigridVariable.solve(_soluteSpecies, _particulateSpecies,
				_boundaryLayer);
		// perform the global reactor mass balances
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].initializeGlobalRateSeries();
		}
	}

	/**
	 * Updates the bulk concentrations of solutes and the discrete
	 * representation of the solids concentration fields and then determines the
	 * concentration fields of all solutes by solving the diffusion/reaction
	 * PDE's using multigrid
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void updateSoluteConcentrations()
			throws MultigridSystemNotSetException {
		// perform the global reactor mass balances
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].updateBulkConcentrationAndRateSeries(_timeStep);
		}
		updateBioDiscreteData();
		solveWithMultigrid();
		// NOTE to solve by relaxation use:
		// MultigridVariable.solveByRelax(
		// _chemicalSpecies,
		// _fixedSpecies,
		// _boundaryLayer);
	}

	/**
	 * Determines the concentration fields of all solutes by solving the
	 * diffusion/reaction PDE's using multigrid
	 * 
	 * @throws MultigridSystemNotSetException
	 * 
	 */
	public void solveWithMultigrid() throws MultigridSystemNotSetException {
		MultigridVariable.solve(_soluteSpecies, _particulateSpecies,
				_boundaryLayer);
	}

	/**
	 * Compute the time step for growth of biomass and global mass balances.
	 * Uses TimeStepConstraint in some operations to compare timesteps and keep
	 * track of which is the constraint setting the time step.
	 */
	private void computeTimeStep() throws BiofilmProblemException {
		TimeStepConstraint presentTimeStepConstraint = new TimeStepConstraint();
		if (!_overrideTimeStep) {
			// ///////determine the variable time step
			// Constraint 1: minimum particle doubling time
			// Has to be the first to be computed!!!
			TimeStepConstraint particleConstraint = biomassContainer
					.getParticleTimeConstraint();

			// Constraint 2: global volume change
			TimeStepConstraint biofilmSizeConstraint = biomassContainer
					.getGlobalSizeConstraint();
			// NOTE comment previous line and uncomment following if no biofilm
			// growth constraint is required
			// float biofilmSizeConstraint = Float.POSITIVE_INFINITY;
			// Constraint 3: the computation of bulk concentration time
			// constraints
			TimeStepConstraint bulkConcentrationConstraint = new TimeStepConstraint();
			for (int i = 0; i < _soluteSpecies.length; i++) {
				bulkConcentrationConstraint = TimeStepConstraint.getMinimum(
						bulkConcentrationConstraint, _soluteSpecies[i]
								.getBulkSpeciesTimeStepConstraint());
			}
			// other timestep limiting criteria may be added here...
			presentTimeStepConstraint = TimeStepConstraint.getMinimum(
					particleConstraint, biofilmSizeConstraint,
					bulkConcentrationConstraint);
			_timeStep = presentTimeStepConstraint.getTimeStep();
			// if time step is infinity it means biofilm is dead!!!
			// break iteration
			if (_timeStep == Float.POSITIVE_INFINITY)
				throw new BiofilmProblemException();
			// If a value for the compulsory time step is defined, determine if
			// compulsory time step should be preformed instead of using the
			// variable time step just determined
			if (_compulsoryTimeStep != 0) {
				// compute the number of compulsory steps already passed
				double nCompulsories = Math.floor(_time / _compulsoryTimeStep);
				// compute the nexte compulsory step
				double nextCompulsory = (nCompulsories + 1)
						* _compulsoryTimeStep;
				// compute the time remaining to next compulsory
				float timeToCompulsory = (float) nextCompulsory - _time;
				if (timeToCompulsory < _timeStep * 1e-3f) {
					nextCompulsory += _compulsoryTimeStep;
					timeToCompulsory += _compulsoryTimeStep;
				}
				// Due to float precision, this operation may cause the
				// simulation to get stuck if timeToCompulsory returns a very
				// small value. Therefore, a check is introduced:
				// _timeStep is used not only if:
				// timeToCompulsory < _timeStep
				// but also if:
				// timeToCompulsory > _timeStep * 1e-3
				// or
				// _writeTimedWriters = false
				if ((timeToCompulsory <= _timeStep)
						| (nextCompulsory == (_time + _timeStep))) {
					// the check (nextCompulsory == (_time + _timeStep))
					// was also included because of precision issues
					_timeStep = timeToCompulsory;
					presentTimeStepConstraint.setTimeStep(_timeStep);
					presentTimeStepConstraint.setName("Compulsory time step");
					// set the flag to write in this iteration
					_writeTimedWriters = true;
				} else {
					// set the flag NOT to write in this iteration
					_writeTimedWriters = false;
				}
			}
		} 
		//ADDITION JONAS TO ACCOMODATE OVERRIDInG TIMESTEPS WITH SPECIFIC TIMES AT WHICH TO OUTPUT AL DATA
			else {
				if (_overrideTimeStepWithSpecificTimesForAllData){
					if(_tOut>0 && _iterationCounter%_tOut==0 ){
						_writeTimedWriters = true;
					} else {
						_writeTimedWriters = false;
					}
					
					
				}	
		}
		// If a maximum time step is defined, truncate at that value
		if (_maximumTimeStep > 0) {
			if (_maximumTimeStep < _timeStep) {
				_timeStep = _maximumTimeStep;
				presentTimeStepConstraint.setTimeStep(_maximumTimeStep);
				presentTimeStepConstraint.setName("Maximum time step");
			}
		}
		// If a minimum time step is defined, lift to at least that value
		if (_minimumTimeStep > 0) {
			_timeStep = Math.max(_timeStep, _minimumTimeStep);
			if (_minimumTimeStep > _timeStep) {
				_timeStep = _minimumTimeStep;
				presentTimeStepConstraint.setTimeStep(_minimumTimeStep);
				presentTimeStepConstraint.setName("Minimum time step");
			}
		}
		// update counters
		incrementIteration();
		_time += _timeStep;
		_timeSeries.add(_time);
		_timeStepConstraints.add(presentTimeStepConstraint);
	}

	/**
	 * Increment the iteration counter by 1
	 */
	public void incrementIteration() {
		_iterationCounter++;
	}

	/**
	 * Perform the growth and division part of the iterative step
	 */
	public void performGrowthAndDivision() throws ModelException {
		// reset the reaction rates for all the reactions
		for (Iterator<Reaction> iter = _reactions.iterator(); iter.hasNext();) {
			Reaction r = iter.next();
			r.resetGlobalReactionRate();
		}
		// compute the timestep for present iteration
		computeTimeStep();
		// Here are the growth steps
		// 0. Shuffle entries in particle list.
		biomassContainer.shuffleParticleList();
		// 1. GROWTH step: only updates masses of all bacteria
		biomassContainer.grow(_timeStep);
		// TODO the following is preliminary code for the attachment of biomass
		// dependent o n the gorwth rate of the bacteria species:
		// ((BiomassParticleContainer) biomassContainer).attach();
		// 2. DEATH step: remove cells with mass below a critical level
		biomassContainer.removeDeadBiomass();
		// 3. DIVISION step: splits bacteria over critial mass
		biomassContainer.divideAndExcreteEPS();
		// 4. DETACHMENT step: performs erosion, sloughing and removes
		// out-of-bounds particles
		// detach();
	}

	/**
	 * Spread biomass after growth and division
	 */
	public void spread() {
		biomassContainer.spread();
	}

	public void spreadByShovingCompletely() {
		biomassContainer.spreadByShoving();
	}

	/**
	 * Perform a single spreading step (for testing or illustration)
	 */
	public void performSpreadingStep() {
		biomassContainer.performSpreadingStep();
	}

	/**
	 * Perform the biomass shrinking by pressure
	 */
	public void performSpreadingByPulling() {
		biomassContainer.performSpreadingByPulling();
	}

	/**
	 * Remove all biomass to be dettached
	 */
	public void detach() throws ModelException {
		biomassContainer.removeDetachedBiomass();
	}

	/**
	 * Get the biomass at present state as collection of Particles
	 * 
	 * @return collection of biomass particles
	 */
	public Collection getBiomassAsBiomassParticleCollection() {
		return biomassContainer.getBiomassAsBiomassParticleCollection();
	}

	/**
	 * Write the concetrations of all solutes to text files in a subdirectory of
	 * directory dirName. The subdirectory 'iteration*' is created, where * is a
	 * variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeSoluteConcentrationsToFile(String dirName)
			throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].finestGridToFile(dir.toString() + "/");
		}
	}

	/**
	 * Write the concetrations of all solids to text files in a subdirectory of
	 * directory dirName. The subdirectory 'iteration*' is created, where * is a
	 * variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeSolidConcentrationsToFile(String dirName) throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].finestGridToFile(dir.toString() + "/");
		}
	}

	/**
	 * Write the values in a multigrid variable to text files in a subdirectory
	 * of directory dirName. The subdirectory 'iteration*' is created, where *
	 * is a variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeMultigridVariableToFile(String dirName,
			MultigridVariable var) throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		var.finestGridToFile(dir.toString() + "/");
	}

	/**
	 * Read the discrete concentrations of all solutes from file
	 * 
	 * @param inputDirectory
	 *            the directory where solute concentrations are
	 * @param inputIterationNumber
	 *            the iteration number to input
	 */
	public void readBiomassMatricesFromDisk(String inputDirectory,
			int inputIterationNumber) {
		String dir = inputDirectory + "/" + SolidsConcentrationWriter.SOLIDDIR
				+ "/iteration" + getFormatedNumber(inputIterationNumber);
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].readDiscreteConcentrationFromDir(dir,
					inputIterationNumber);
		}
	}

	public void readBulkConcentrationsFromDisk(String inputDirectory,
			int inputIterationNumber) {
		String resultFile = inputDirectory + "/"
				+ SimulationResultsWriter.FILENAME;
		SimulationResultParser p = new SimulationResultParser(resultFile);
		for (int i = 0; i < _soluteSpecies.length; i++) {
			float v = p.getValueForBulkConcentration(_soluteSpecies[i]);
			BulkConcentration bc = new ConstantBulkConcentration(v);
			_soluteSpecies[i].setBulkConcentration(bc);
		}

	}

	/**
	 * Writes the detachment level set matrix to file, before detachment is
	 * carried out
	 * 
	 * @param dirName
	 * @throws Exception
	 */
	public void writeDetachmentLevelSet(String dirName) throws Exception {
		biomassContainer.writeDetachmentLevelSet(dirName + "/iteration"
				+ getFormatedIterationNumber());
	}

	/**
	 * Write the growing cells of a biofilm in boolean matrix to file
	 * 
	 * @param dirName
	 * @throws IOException
	 *             if write to file fails
	 */
	public void writeBiofilmFront(String dirName) throws IOException {
		biomassContainer.writeBiofilmFront(dirName + "/iteration"
				+ getFormatedIterationNumber());
	}

	/**
	 * Writes the positions of bacteria as particles
	 * 
	 * @param dirName
	 *            the directory to write text file to
	 * @throws IOException
	 *             in case of failure to write file
	 */
	public void writeParticlePositionsToFile(String dirName) throws IOException {
		biomassContainer.bacteriaToFile(dirName + "/iteration"
				+ getFormatedIterationNumber() + ".txt");
	}

	/**
	 * @return the present time
	 */
	public float getTime() {
		return _time;
	}

	/**
	 * @return the present iteration number
	 */
	public int getIterationCounter() {
		return _iterationCounter;
	}

	/**
	 * @return the present iteration number formated to a String
	 */
	public String getFormatedIterationNumber() {
		return getFormatedNumber(_iterationCounter);
	}

	/**
	 * @param n
	 *            number to return in formatted string
	 * @return a number formated to a String
	 */
	public String getFormatedNumber(int n) {
		NumberFormat numberFormat = new DecimalFormat("########");
		numberFormat.setMinimumIntegerDigits(8);
		return numberFormat.format(n);
	}

	/**
	 * Get the current number of particles of a given species
	 * 
	 * @param s
	 * @return number of particles of a given species
	 */
	public int getNumberOfParticles(BiomassSpecies s) {
		return biomassContainer.getNumberOfParticles(s);
	}

	/**
	 * Return a tab separated string with some parameters concerning the present
	 * iteration: iteration number, time
	 * 
	 * @return string with some parameters
	 */
	public String getIterationParameters() {
		long t = (System.currentTimeMillis() - _realInitialTime);
		// The next line is to include real time of simulations
		String s = _iterationCounter + "\t" + _time + "\t" + t;
		return s;
	}

	/**
	 * @return the header for a table with iteration parameters
	 */
	public String getIterationParametersHeader() {
		// The next line is to include real time of simulations
		return "Iteration\tTime [h]\tCurrent real time [ms]";
	}

	/**
	 * @return Returns the dimensionality of the system (2 for 2D, 3 for 3D)
	 */
	public int getDimensionality() {
		return _dimensionality;
	}

	/**
	 * @return Returns the carrierArea.
	 */
	public float getComputationalVolumeMultiplier() {
		return _carrierArea / getComputationalVolumeCarrierArea();
	}

	/**
	 * Set the ratio between carrier area in the reactor and the computational
	 * volume carrier area
	 * 
	 * @param m
	 *            The computational volume multiplier.
	 */
	public void setCarrierArea(float m) {
		this._carrierArea = m;
	}

	/**
	 * @return Returns the residenceTime.
	 */
	public final float getResidenceTime() {
		return _residenceTime;
	}

	/**
	 * @param residenceTime
	 *            The residenceTime to set.
	 */
	public void setResidenceTime(float residenceTime) {
		this._residenceTime = residenceTime;
	}

	/**
	 * @return Returns the _reactorVolume.
	 */
	public float getReactorVolume() {
		return _reactorVolume;
	}

	/**
	 * @param volume
	 *            The _reactorVolume to set.
	 */
	public void setReactorVolume(float volume) {
		_reactorVolume = volume;
	}

	/**
	 * @return Returns the time series of simulation.
	 */
	public UnlimitedFloatArray getTimeSeries() {
		return _timeSeries;
	}

	/**
	 * @return Returns an unmodifiable collection with solute species.
	 */
	public Collection<SoluteSpecies> getSoluteSpecies() {
		return Collections
				.unmodifiableCollection(Arrays.asList(_soluteSpecies));
	}

	/**
	 * @return the number of solute species in the system
	 */
	public int getNumberOfSolutes() {
		return _soluteSpecies.length;
	}

	/**
	 * Solute species at index i in the _soluteSpecies array
	 * 
	 * @param i
	 * @return the solute species number i
	 */
	public SoluteSpecies getSoluteSpecies(int i) {
		return _soluteSpecies[i];
	}

	/**
	 * @return Returns an unmodifiable collection with particulate species.
	 */
	public Collection<ParticulateSpecies> getParticulateSpecies() {
		return Collections.unmodifiableCollection(Arrays
				.asList(_particulateSpecies));
	}

	/**
	 * @return Returns an unmodifiable collection with biomass species.
	 */
	public Collection<BiomassSpecies> getBiomassSpecies() {
		return Collections.unmodifiableCollection(Arrays
				.asList(_biomassSpecies));
	}

	/**
	 * @return Returns the _detachedBiomass.
	 */
	public DetachedBiomassContainer detachedBiomassContainer() {
		return _model._detachedBiomass;
	}

	/**
	 * @return Returns the biomassContainer.
	 */
	public BiomassContainer getBiomassContainer() {
		return biomassContainer;
	}

	/**
	 * @return true if (_time >= _finishIterationTime) or _finishSimulation
	 */
	public boolean endSimulation() {
		
		// otherwise, check other criteria
		//Jonas
		if (!(_nutrientForConsumptionStopCriterion == null)){
			
		return (_time >= _finishIterationTime) | _finishSimulation
				| (getCurrentBiovolume() >= _maxBiovolume)
				| (getCurrentBiofilmHeight() > _finalBiofilmHeight
						// Jonas
						| stopDueToConsumption(_nutrientForConsumptionStopCriterion));
		}
		//JONAS. TO remove JONAS action take the following out of th eelse loop and replace everything between //Jonas and 
		//end JONAS with it
		else {
			return (_time >= _finishIterationTime) | _finishSimulation
					| (getCurrentBiovolume() >= _maxBiovolume)
					| (getCurrentBiofilmHeight() > _finalBiofilmHeight);
			
			}	
		//end Jonas
		
	}

	/**
	 * Update values of maximum specific growth rates for each fixed species
	 */
	public void updateValuesOfCurrentMaximumSpecificGrowthRates() {
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].computeMaximumSpecificGrowthRate();
		}
	}

	/**
	 * @return the area represented in the computational volume [um^2}
	 */
	public float getComputationalVolumeCarrierArea() {
		return (_dimensionality == 2 ? systemSize.y * _systemGridSide
				: systemSize.y * systemSize.z);
	}

	/**
	 * @return Returns the depth of the system for 2D computations.
	 */
	public float get2DSystem3rdDimension() {
		if (_dimensionality == 3)
			throw new ModelRuntimeException(
					"Trying to get System depth in a 3D system");
		return _systemGridSide;
	}

	/**
	 * @return the side of a grid element
	 */
	public float getGridElementSide() {
		return _systemGridSide;
	}

	/**
	 * Set a time criterium for finishing simulations. Default is
	 * Float.POSITIVE_INFINITY, meaning that simulations must be stopped
	 * manually
	 * 
	 * @param iterationTime
	 *            The _finishIterationTime to set.
	 */
	public void setFinishIterationTime(float iterationTime) {
		_finishIterationTime = iterationTime;
	}

	/**
	 * Set the iteration time step to a fixed value, overriding the models
	 * variable time step feature
	 * 
	 * @param t
	 *            time step
	 */
	public void overrideTimeStep(float t) {
		_timeStep = t;
		_overrideTimeStep = true;
		_writeTimedWriters = true;
	}
	/**
	 * Set the iteration time step to a fixed value, overriding the models
	 * variable time step feature
	 * @addition by Jonas: sometimes you may want to override the timestep but output data less frequently. 
	 * this is the second input here.
	 * @param t
	 *            time step
	 * @param tOut
	 * 		      every tOut timesteps: output all data
	 */
	public void overrideTimeStepWithSpecificTimesForAllData(float t, int tOut) {
		_overrideTimeStepWithSpecificTimesForAllData = true;
		_timeStep = t;
		_overrideTimeStep = true;
		_writeTimedWriters = false;
		_tOut = tOut;
	}
	/**
	 * Turn off the overriding of time step
	 */
	public void setTimeStepBackToManual() {
		_timeStep = 0;
		_overrideTimeStep = false;
		_writeTimedWriters = true;
	}

	/**
	 * Set a value for the compulsory time step. This means that at each time
	 * interval defined by attribute _compulsoryTimeStep an iteration will be
	 * performed, inspite of the variable time step. This feature may be used so
	 * that output is presented at regular time intervals
	 * 
	 * @param f
	 */
	public void setCompulsoryTimeStep(float f) {
		_compulsoryTimeStep = f;
		_overrideTimeStep = false;
	}

	/**
	 * Sets a maximum time step, writing outputs at each iteration independently
	 * of the time step value
	 * 
	 * @param f
	 */
	public void setMaximumTimeStep(float f) {
		_maximumTimeStep = f;
	}

	/**
	 * @return Returns true if timed state writers should be written at this
	 *         iteration.
	 */
	public boolean writeTimedWriters() {
		return _writeTimedWriters;
	}

	/**
	 * set the simulation to finish
	 */
	public void setFinishSimulation() {
		_finishSimulation = true;
	}

	/**
	 * @return Returns the _boundaryLayer.
	 */
	public BoundaryLayer getBoundaryLayer() {
		return _boundaryLayer;
	}

	/**
	 * @param t
	 *            The value for the minimum time step in the model simulation.
	 */
	public void setMinimumTimeStep(float t) {
		_minimumTimeStep = t;
	}

	/**
	 * @return the last time constraint or a default time contraint
	 */
	public TimeStepConstraint getLastTimeConstraint() {
		int l = _timeStepConstraints.size();
		if (l > 0) {
			TimeStepConstraint c = (_timeStepConstraints.get(l - 1));
			return c;
		}
		// if iterations have not started yet, return a new TimeConstraint
		return new TimeStepConstraint(0, "Initial time point");
	}

	/**
	 * @return the Singleton instance of model for dessirialization
	 */
	protected Object readResolve() {
		_model = this;
		return _model;
	}

	private BiomassBooleanMatrix _biomassBooleanMatrix;

	/**
	 * Returns a boolean matrix with the locations of the biomass
	 * 
	 * @return the present instance of BiomassBooleanMatrix
	 */
	public BiomassBooleanMatrix getBiomassBooleanMatrix() {
		if (_biomassBooleanMatrix == null)
			_biomassBooleanMatrix = new BiomassBooleanMatrix();
		_biomassBooleanMatrix.upateBiomassConcentrations(_particulateSpecies);
		return _biomassBooleanMatrix;
	}

	/**
	 * Set a maximum biovolume as a criterium to end simulations
	 * 
	 * @param biovolume
	 */
	public void setMaxBiovolume(float biovolume) {
		_maxBiovolume = biovolume;
	}

	/**
	 * Set a maximum runlength as a criterium to end simulations. Only works for
	 * granule
	 * 
	 * @param r
	 *            maximum run length
	 */
	public void setMaxRunLength(float r) {
		_maxRunLength = r;
	}


	/**
	 * Set maximum biofilm height as criterium to end simulaiton.
	 * 
	 * @param maxHeight
	 */
	public void setMaxHeight(float maxHeight) {
		_finalBiofilmHeight = maxHeight;
	}
	
	//Jonas
	/**
	 * Set maximum biofilm height as criterium to end simulaiton.
	 * 
	 * @param maxHeight
	 */
	public void setNutrientForConsumptionStopCriterion(String s) {
		for (int i = 0; i < _soluteSpecies.length; i++) {
			if(_soluteSpecies[i].getName()==s){		
					_nutrientForConsumptionStopCriterion = _soluteSpecies[i];
					}
		}

	}
	
	public void setMaxConsumption(float c) {
		_maxConsumption = c;
	}
	public boolean stopDueToConsumption(SoluteSpecies s){
		_cummulativeConsumptionOfNutrient += s.computeGlobalRateFromMassBalance();
		if (_cummulativeConsumptionOfNutrient<=_maxConsumption){
			return true;
		}
		return false;
		
	}
	//end JOnas
}