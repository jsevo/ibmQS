/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * A class that derives from reaction and that implements and inlfux or outflux
 * of a given solute species, the one set as attribute _species. This class is
 * used to model fluxes into "birds-eye-view" models.
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class Flux extends Reaction implements Serializable {
	private SoluteSpecies _species;

	private float _bulkConcentration;

	/**
	 * creates a new reaction with one or more factors
	 * 
	 * @param n
	 *            name of reaction
	 * @param bs the solute species that is flowing in or out
	 * @param bulkConcentration
	 */
	public Flux(String n, SoluteSpecies bs, float bulkConcentration) {
		super(n, bs, 0, 0);
		_species = bs;
		_bulkConcentration = bulkConcentration;
		// add the reaction thus created to the Model._reactions
		Model.model().addReaction(this);
	}

	/**
	 * Returns the rate at the present location
	 * 
	 * @return flux rate [g/um^3/h]
	 */
	public float getRate() {
		return (_bulkConcentration - _species.getValue());
	}

	/**
	 * Pre-compute the mass-based growth rate and store it in atribute
	 * _presentReactionRate. Also, update the value of _globalReactionRate
	 * 
	 * @param mCatalyst
	 *            mass of catalyst species
	 */
	public void computeMassGrowthRateAndAddToGlobal(BiomassSpecies.Composition c) {
		// overrides original but does nothing, since the influx will never
		// be used for kinetics of particulate species
	}

	/**
	 * Just returns 0. Never used in particulate species rates
	 * 
	 * @return the mass growth rate [g/h]
	 */
	public float getPreComputedMassGrowthRate() {
		return 0;
	}

	/**
	 * Just returns 0. Never used in particulate species rates
	 * 
	 * @return specific rate of reaction [h^-1]
	 */
	public float getSpecificRateFactor() {
		return 0;
	}

	/**
	 * Just returns 0. Never used in particulate species rates
	 * 
	 * @param c
	 *            composition of biomass particle to compute specific rate
	 * @return specific rate of reaction [h^-1]
	 */
	public float getSpecificRateFactor(BiomassSpecies.Composition c) {
		return 0;
	}

	/**
	 * Just returns 0. Never used in particulate species rates
	 * 
	 * 
	 * @return maximum specific rate in the system [h^-1]
	 */
	public float getSpecificRateMaximum() {
		return 0;
	}

	/**
	 * Returns the derivative of the flux rate at the present
	 * location
	 * 
	 * @param c
	 *            chemical species to derivate rate to
	 * @return derivative of rate of reaction [g/um^3/h]
	 */
	public float getRateDerivative(SoluteSpecies c) {
		if (c == _species) {
			return -_species.getValue();
		}
		// if the species is not the catalyst
		return 0;
	}

	/**
	 * Update the array with values of the rate and rate derivative
	 * 
	 * @param c
	 *            the solute species for which this is being computed
	 * @param rDr
	 *            [rate, rateDerivative]
	 */
	public void updateValuesForRateAndRateDerivative(SoluteSpecies c,
			float[] rDr) {
		rDr[0] = getRate();
		rDr[1] = getRateDerivative(c);
	}

}