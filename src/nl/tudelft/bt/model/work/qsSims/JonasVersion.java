	package nl.tudelft.bt.model.work.qsSims;

import java.awt.Color;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.BiomassVizualizer;
import nl.tudelft.bt.model.apps.output.BiofilmMaximumThicknessSeries;
import nl.tudelft.bt.model.apps.output.FixedTotalBiomassSeries;
import nl.tudelft.bt.model.apps.output.ParticlePositionWriter;
import nl.tudelft.bt.model.apps.output.ProducedBiomassSeries;
import nl.tudelft.bt.model.apps.output.SimulationResultsWriter;
import nl.tudelft.bt.model.apps.output.SolidsConcentrationWriter;
import nl.tudelft.bt.model.apps.output.SoluteConcentrationWriter;
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
 * 
 * 
 * @author Jonas Schluter (jonas.schluter@dtc.ox.ac.uk)
 */
public class JonasVersion extends ModelHandler {
	protected static int geometry = 2;
	// ******************** Input Parameters
	// *************************************
	protected static String outputDirectory = "/host/JonasWork/simulations/28022011/security";
	// SPECIFICITY PARAMETERS
	protected static float uMax =1f;
	/*
	 * Two substrates multiplicator: Because I have two substrats between which switching occurs muMax'es for the substrates are included in 
	 * the calculation of the "process factor".  Multiplying with another muMax would just scale up the calculated rate, thus == 1 is suggested. 
	 */
	// DIFFUSION COEFFICIENTS [um2/h]
	private static float Dn = 4e4f; 
    // *****************************************************************************
	// ********************** BULK CONCENTRATIONS
	// *********************************
	protected static float Nbulk = .153f; // nutrien bulk concentration (zero as it comes from the solid medium)
	// **********************specific masses
	// *********************
	protected static float specificMassBa = 150f; // Joao
	// ********************** Finish conditions / other system parameters
	static float outputEvery = 0.05f;
	static float simulationFinishTime = 80.0f; // This will be limiting
	// ********************** Yields and Rates *********************
	// Yield coefficients (positive)
	private static float YA = 0.5f; // Sara, Jonas: [gX/gS] 
	private static float KSGrowthBa = 3.5e-5f; // Sara
	// **************Computation parameters**************************************
	protected static float systemSize = 200f; // Joao, Sara [um]
	protected static float relativeMaximumRadius = 1f / systemSize;
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f; // Joao
	protected static float relativeBoundaryLayer = (float) 25 / systemSize; // Joao
	protected static int gridSide = 65; // Joao, multigrid grid side
	protected static float kShov = 1.0f; // Joao, shoving parameter
	protected static float rdetach = 1f;



	// *********************** Initial Particle numbers etc.
	// ***********************
	// initial number of particles in the system (inoculum)
	protected static int iA;//initial particle number A
	protected static int iB;
	private static BiomassSpecies biomassA; // 
	private static BiomassSpecies biomassB; // 
	/**
	 * Define the bacteria species, the chemical species and the processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {

		// 1. CREATE THE SOLUTES
		// ###################################################
		// create the solute species
		SoluteSpecies N = new SoluteSpecies("N", Dn);
		SoluteSpecies pG = new SoluteSpecies("pG", Dn);
		SoluteSpecies iN = new SoluteSpecies("iN", Dn);
		// set bulk concentrations
		N.setBulkConcentration(new ConstantBulkConcentration(Nbulk));	
		pG.setBulkConcentration(new ConstantBulkConcentration(0f));
		iN.setBulkConcentration(new ConstantBulkConcentration(1f));
		// 2. CREATE THE PARTICULATE SPECIES
		// #######################################
		ParticulateSpecies activeA = new ParticulateSpecies("A", specificMassBa, Color.green);
		ParticulateSpecies activeB = new ParticulateSpecies("B", specificMassBa, Color.blue);
		// array of fixed species that constitute biomass species
		ParticulateSpecies[] spA = { activeA };

		float[] fractionalVolumeCompositionA = { 1.0f };

		// 3. CREATE THE BIOMASS
		// SPECIES############################################
		//good guys
		biomassA = new BiomassSpecies("speciesA",
				spA, fractionalVolumeCompositionA);
		biomassA.setActiveMass(activeA);
		//biomassA.setInducibleColor(N, KSGrowthBa/10, Color.RED, Color.green);
		
		// 4. CREATE THE REACTION FACTORS, MONOD AND INHIBITION
		// COEFFICIENTS########
		ProcessFactor GrowthN = new Saturation(N,KSGrowthBa);//, (uMaxGrowthGoodGuys>uMaxfucoseUtilizationGoodGuys)); // bacterial growth on normal Nutrients (from above).
		ProcessFactor GrowthStep = new Step(iN,.9f);
		// 5. CREATE THE REACTIONS
		// #################################################
		// Normal growth, nutrients from 'above'
		Reaction growthAN = new Reaction("growthAN",
				activeA, uMax, 1);
		growthAN.addFactor(GrowthN);
		
		Reaction pGoodProd = new Reaction("pGoodProd", activeA, uMax,1);
		pGoodProd.addFactor(GrowthN);
		//pGoodProd.addFactor(GrowthStep);
		
		// 6. ASSIGN REACTION TO THE SPECIES THROUGH REACTION
		// STOICHIOMETRIES#######
		// 
		// active mass
		NetReaction rsA = new NetReaction(2);
		rsA.addReaction(growthAN, YA);
		rsA.addReaction(pGoodProd, 0);
		activeA.setProcesses(rsA);



		// assign reaction stoichiometry to the solutes
		// solutes
		NetReaction rsN = new NetReaction(1);
		rsN.addReaction(growthAN, -1/YA);
		N.setProcesses(rsN);

		NetReaction rsP = new NetReaction(1);
		rsP.addReaction(pGoodProd, 1);
		pG.setProcesses(rsP);

		NetReaction rsIn = new NetReaction(0);
		iN.setProcesses(rsIn);
		// 7. ADD THE SOLUTE SPECIES AND THE BIOMASS SPECIES (which contain the
		// particulate species) TO THE SYSTEM
		// ######################################
		addBiomassSpecies(biomassA);
		addSoluteSpecies(N);
		addSoluteSpecies(pG);
		addSoluteSpecies(iN);
	}

	// #############################inoculate######################################
	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	protected void inoculate() {
/*
		int totalCellsToPlace = iA;
		for (int i = 0; i < totalCellsToPlace; i++) {
			// set the center to a random position at the surface
			float y = Model.model().getRandom();
			float z = Model.model().getRandom();
				placeBiomass(biomassA,mucusHeight,z*systemSize,0);
		}
	}

*/
		int[] ncells = {iA};
		inoculateRandomly(ncells);
	}


	/*
	 * INITIALIZE THE DETACHMENT
	 * FUNCTION########################################### (non-Javadoc)
	*/
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new Height2MassDetachment(0f);
		setDetachmentHandler(df);
	}
 /*
	
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new OnlyOutOfBoundsDetachment(1);
		setDetachmentHandler(df);
	}
	*/
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
							+ "1: output directory (CAUTION!!! directory will be erased \n"
							+ "2: random seed"
							+ "3: iA"
							+ "4: iB");
		}

		outputDirectory = args[0];
		iA = (int) Float.parseFloat(args[2]);
		iB = (int) Float.parseFloat(args[3]);
		boolean allData = true;
		boolean visualsOn = true;
		int seed = Integer.parseInt(args[1]);
		
		Model.model().setSeed(seed);
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new JonasVersion();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		if (visualsOn){
		 app = new BiomassVizualizer(app);
		 //app = new VizualModelControler(app);
		}
		 try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			if (allData){
			//app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			}
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);

			app.addStateWriter(spw);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			/*for (Iterator i = Model.model().detachedBiomassContainer()
					.getCustomizedSeries().iterator(); i.hasNext();) {
				DetachedSeries s = (DetachedSeries) i.next();
				spw.addSeries(s);
			}*/
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
			// start iterating cycle

			// ***** FINISH SIMULATION CONDITIONS *****

			Model.model().setFinishIterationTime(simulationFinishTime);


		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// grow once to make sure cells are colored when firs writing
			//Model.model().overrideTimeStep(0.001f);
			//app.writeState();
			//app.performFullIteration();
			//app.writeState();

			// start iterating cycle
			//Model.model().overrideTimeStep(outputEvery);
			app.startIterating();
			System.out.println("Simulation finished.");
			System.out.println(Model.model().getTimeStep());
			System.exit(-1);
		} catch (Exception e1) {
			e1.printStackTrace();
		}

	}
}