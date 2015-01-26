package nl.tudelft.bt.model.work.qsSims;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.BiomassVizualizer;
import nl.tudelft.bt.model.apps.output.BiovolumeSeries;
import nl.tudelft.bt.model.apps.output.FixedTotalBiomassSeries;
import nl.tudelft.bt.model.apps.output.ParticlePositionWriter;
import nl.tudelft.bt.model.apps.output.PovRayWriter;
import nl.tudelft.bt.model.apps.output.ProducedBiomassSeries;

import nl.tudelft.bt.model.apps.output.SimulationResultsWriter;
import nl.tudelft.bt.model.apps.output.SolidsConcentrationWriter;
import nl.tudelft.bt.model.apps.output.SoluteConcentrationWriter;
import nl.tudelft.bt.model.apps.output.TimedStateWriterDecorator;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.reaction.NetReaction;
import nl.tudelft.bt.model.reaction.ProcessFactor;
import nl.tudelft.bt.model.reaction.Reaction;
import nl.tudelft.bt.model.reaction.Saturation;
import nl.tudelft.bt.model.reaction.Step;
/**
 * Simulates growth of two strains with differential public good production
 * in a vertical slice biofilm cross-section
 * 
 * @author Jonas Schluter (jonas.schluter+github@gmail.com)
 */
public class new4_wideSystem_inducerProductionVarINI_matchQN_viaRE extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory;


	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	
	//protected static float maxHeight = 120f; // [um]
	protected static float simulationFinishTime = 80; // [h]
	
	protected static int totalNumber = 0;
	protected static float ratioOfFirstSpecies = 0.5f;


	// 1. Reaction parameters
	// *********************************************************************************

	// max growth rate on substrate alone
	// private static float uMax = 1.5f;
	
	private static float uMax = 1f;

	private static float Ks = 3.5e-5f;
	
	// Effect of public goods upon growth rate
	private static float benefit = 3f; // 1.1 //[1/T
	
	// Threshold public good concentration before benefit is received
	private static float goodThreshold = 4e-3f; // 4e-3f
	
	// Threshold inducer concentration before public good production is started
	private static float inducerThreshold = 1f;
	
	// Rate of goods production
	private static float goodRate1 = 0.8123f; // [1/T]
	
	// Rate of inducer production
	private static float inducerRate1 = 1f; // [1/T]

	private static float costScale = .3f; // 1e-5f; // scales the cost of the

	// public good

	// Yield of biomass on substrate
	private static float Yxs = 0.5f; // [gCOD-PHB/gCOD-S]

	// //////////////////////////////////////////////////////////////////////////////////
	// Solute and Biomass Parameters

	// Public Good
	protected static float goodBulkConcentration = 0f; // [gGood/L]
	protected static float goodDiffusivity = 3e5f;   // [um2/h]

	// Substrate (Bulk concentration set by input)
	protected static float substrateBulkConcentration = 0.125f;
	private static float substrateDiffusivity = 4e4f; // [um2/h]
	
	// Inducer
	protected static float inducerBulkConcentration = 0f; // [gGood/L]
	protected static float inducerDiffusivity = 3e5f;   // [um2/h]

	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// //////////////////////////////////////////////////////////////////////////////////
	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 400; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 1f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = (float) 25 / systemSize;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber1 = 0;

	protected static int initialParticleNumber2 = 0;

	// outpute (write results to file) every:
	protected static float outputEvery =  .1f; //[h]
	private static float finishIterationTime = 48f;
	//protected static float maxBiovolume = 40000f; // [uL]


	// /END OF PARAMETERS

	/**
	 * Define bacterial species, solute species, and reaction processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// 2a. Create the solute species (public good only)
		// *******************************************************************************
		// substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		// public good
		SoluteSpecies pGood = new SoluteSpecies("pGood", goodDiffusivity);
		pGood.setBulkConcentration(new ConstantBulkConcentration(
				goodBulkConcentration));
		
		// inducer
		SoluteSpecies inducer = new SoluteSpecies("inducer", inducerDiffusivity);
		inducer.setBulkConcentration(new ConstantBulkConcentration(
				inducerBulkConcentration));
		
		//public ItermittentBulkConcentration(float vfeast, float vfamine,
		//float initalFeastTime, float feastTime, float famineTime)

		// *******************************************************************************
		// Strain One, public good secreting
		// Active mass of strain 1, color to be overridden
		ParticulateSpecies activeOne = new ParticulateSpecies("activeOne",
				specificMassBiomass, Color.blue);
		
		ParticulateSpecies[] spOne = { activeOne };

		float[] fractionalCompositionOne = { 1.0f };

		BiomassSpecies speciesOne = new BiomassSpecies("speciesOne", spOne,
				fractionalCompositionOne);
		speciesOne.setActiveMass(activeOne);

		// ///////////////////////////////////////////////////////////////////////
		// Strain two, exploitative
		// Active mass of strain 2
		ParticulateSpecies activeTwo = new ParticulateSpecies("activeTwo",
				specificMassBiomass, Color.red);
		
		ParticulateSpecies[] spTwo = { activeTwo }; 
		
		float[] fractionalCompositionTwo = { 1.0f };

		BiomassSpecies speciesTwo = new BiomassSpecies("speciesTwo", spTwo,
				fractionalCompositionTwo);
		speciesTwo.setActiveMass(activeTwo);

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// *****************************************************************************

		ProcessFactor substrateUtil = new Saturation(substrate, Ks);
		ProcessFactor goodUtil = new StepReturnSubstrateAndGood(substrate, pGood, Ks,
				goodThreshold);
		ProcessFactor inducerStep = new Step(inducer, inducerThreshold);

		// 5. Create growth and goods production reactions for each strain
		// *****************************************************************************
		// Strain 1
		// Growth
		Reaction growthOne;
		growthOne = new Reaction("growthOne", activeOne, uMax, 1);
		growthOne.addFactor(substrateUtil);
	

		// Public Good Production
		Reaction goodProductionOne;
		goodProductionOne = new Reaction("goodProductionOne", activeOne,
				goodRate1, 2);
		goodProductionOne.addFactor(substrateUtil);
		goodProductionOne.addFactor(inducerStep);
	

		// Public good use
		Reaction goodUseOne;
		goodUseOne = new Reaction("goodUseOne", activeOne, benefit*uMax, 1);
		goodUseOne.addFactor(goodUtil);
		
		
		// Inducer production 
		Reaction inducerProductionOne;
		inducerProductionOne = new Reaction("inducerProductionOne", activeOne,
				inducerRate1, 0);
		
		
		
		// ////////////////////////////////////////////////////////////////////////////

		// Strain 2
		// Growth - alternative forms for nutrient graidents ON or OFF
		Reaction growthTwo;
			growthTwo = new Reaction("growthTwo", activeTwo, uMax, 1);
			growthTwo.addFactor(substrateUtil);


		// Public Good Use
		Reaction goodUseTwo;
		goodUseTwo = new Reaction("goodUseTwo", activeTwo, benefit*uMax, 1);
		goodUseTwo.addFactor(goodUtil);
		

		// 6. Assign reactions to the species through ReactionStoichiometries
		// ******************************************************************************
		// Strain 1
		NetReaction rsActiveOne = new NetReaction(3);
		rsActiveOne.addReaction(growthOne, 1);
		rsActiveOne.addReaction(goodUseOne, 1);
		rsActiveOne.addReaction(goodProductionOne, -costScale*uMax);
		activeOne.setProcesses(rsActiveOne);

		// Strain 2
		NetReaction rsActiveTwo = new NetReaction(2);
		rsActiveTwo.addReaction(growthTwo, 1);
		rsActiveTwo.addReaction(goodUseTwo, 1);
		activeTwo.setProcesses(rsActiveTwo);

		// assign reaction stoichiometry to the solute(s)
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growthOne, -1 / Yxs);
		rsSubstrate.addReaction(growthTwo, -1 / Yxs);
		substrate.setProcesses(rsSubstrate);

		// public good (G)
		NetReaction rsPgood = new NetReaction(1);
		rsPgood.addReaction(goodProductionOne, 1);
		pGood.setProcesses(rsPgood);
		
		// inducer (I)
		NetReaction rsInducer= new NetReaction(1);
		rsInducer.addReaction(inducerProductionOne, 1);
		inducer.setProcesses(rsInducer);


		// 7. add solute species and biomass species to the system
		addBiomassSpecies(speciesOne);
		addBiomassSpecies(speciesTwo);

		addSoluteSpecies(substrate);
		addSoluteSpecies(pGood);
		addSoluteSpecies(inducer);

	}

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialParticleNumber1, initialParticleNumber2 };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(kdetach);
		setDetachmentHandler(df);
	}

	/**
	 * @param app
	 */
	protected static void setSystemParametersAndInitializeSystemSpace(
			ApplicationComponent app) {
		// create the space
		app.setSystemSpaceParameters(geometry, systemSize,
				relativeMaximumRadius, relativeMinimumRadius,
				relativeBoundaryLayer, gridSide, kShov);
		// initialize
		app.initializeSystemSpace();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line
		if (args.length < 4) {
			throw new RuntimeException(
					"input arguments missing: \n"
							+ "1: output directory (CAUTION!!! directory will be erased) \n"
							+ "2: seed for random number generator \n"
							+ "3: flag for running with graphics (1 on, 0 off) \n"
							+ "4: initial ratio of first species"
							+ "5: inducer threshold"
							+ "6: initial cell number"
							+ "7: substrateBulkConcentration"
							+ "8: finish iteration time (h)");
		}
				
		// parse inputs
		outputDirectory = args[0];
		int seed = Integer.parseInt(args[1]);
		Model.model().setSeed(seed);
		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
		ratioOfFirstSpecies = Float.parseFloat(args[3]);
		inducerThreshold = Float.parseFloat(args[4]);
		totalNumber = Integer.parseInt(args[5]);
		substrateBulkConcentration = Float.parseFloat(args[6]);
		finishIterationTime = Float.parseFloat(args[7]);

		goodRate1 = 1;
	
		initialParticleNumber1 = (int) (totalNumber * ratioOfFirstSpecies);
		initialParticleNumber2 = (int) (totalNumber * (1f - ratioOfFirstSpecies));
		
		

		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new new4_wideSystem_inducerProductionVarINI_matchQN_viaRE();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VariableSeries biovolume = new BiovolumeSeries();
		

		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// bulk concentrations
			//app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			//app = new VizualModelControler(app);
		}
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
			
			
			// additional output documentation
			app.addStateWriter(new TimedStateWriterDecorator(
					new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			
			/*
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			*/
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// add the time constraints writer
			//app.addStateWriter(new TimeConstraintsWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// add bulk concentrations of all solutes as variable series
			for (Iterator i = Model.model().getSoluteSpecies().iterator(); i
					.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) i.next();
				spw.addSeries(s.getBulkConcentrationSeries());
				spw.addSeries(s.getRateTimeSeries());
			}
			// add particulate global masses to write
			for (Iterator i = Model.model().getParticulateSpecies().iterator(); i
					.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) i.next();
				spw.addSeries(s.getTotalMassSeries());
			}
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			//Model.model().setCompulsoryTimeStep(outputEvery);
			//Model.model().setMaxBiovolume(maxBiovolume);
			//simulationFinishTime = 500; // [h]
			// set a maximum size for ball
			//Model.model().setMaxBiovolume(2e4f);
			//Model.model().setMaxHeight(maxHeight);
			 Model.model().setFinishIterationTime(finishIterationTime);
			 Model.model().setMaxBiovolume(8e4f); //
			 Model.model().overrideTimeStepWithSpecificTimesForAllData(0.15f,400);
			 Model.model().setNutrientForConsumptionStopCriterion("substrate");
			 Model.model().setMaxConsumption(-10e7f);
			 

			// start the iteration
			app.writeState(); // write iteration 0
			app.startIterating();
			//app.forceWriteState();
		} catch (Exception e1) {
			//try {
			//	app.forceWriteState();
			//} catch (IOException e2) {
			//	System.err.println("Error serializing state:");
			//	System.err.println("");
			//	e2.printStackTrace();
			//}
			System.err.println("");
			System.err.println("Program failed due to :");
			System.err.println("");
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}
