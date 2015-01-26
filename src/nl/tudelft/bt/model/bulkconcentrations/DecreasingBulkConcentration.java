/*
 * Created on 10-may-2013 by Sara Mitri
 */
package nl.tudelft.bt.model.bulkconcentrations;

import nl.tudelft.bt.model.Model;

/**
 * Implements bulk concentration of a solute that is reduced over time. 
 * At the begining of each cycle the bulk concentration of the solute is set 
 * to the value of inputConcentration. At each timestep, the concentration of 
 * the solute decreases.
 * 
 * @author Sara Mitri based on SbrBulkConcentration
 */
@SuppressWarnings("serial")
public class DecreasingBulkConcentration extends BulkConcentration {
	private float _inputConcentration;

	private float _maxFractionToDecrease;
	
	private static float DEFAULTMAXFRACTION = 0.95f;

	private float _precisionInConcentration;

	private static float DEFAULTPRECISION = 0.00001f;

	/**
	 * SBR cycles with duration cycleDuration
	 * 
	 * @param inputConcentration
	 *            The starting solute concentration 
	 * @param fractionToDecrease
	 *            the fraction of the liquid in the reaction decreased
	 *            at each timestep
	 * @param precisionInConcentration
	 *            the precision in the concentration value (below this value,
	 *            concentration is assumed 0)
	 */
	public DecreasingBulkConcentration(float inputConcentration) {
		this(inputConcentration, 
				DEFAULTPRECISION,
				DEFAULTMAXFRACTION);
	}

	/**
	 * @param inputConcentration
	 *            The starting solute concentration 
	 * @param fractionToDecrease
	 *            the fraction of the liquid in the reaction decreased
	 *            at each timestep
	 * @param precisionInConcentration
	 *            the precision in the concentration value (below this value,
	 *            concentration is assumed 0)
	 */
	public DecreasingBulkConcentration(float inputConcentration, 
			float precisionInConcentration, 
			float maxFractionToDecrease) {
		super();
		// assign
		_inputConcentration = inputConcentration;
		_maxFractionToDecrease = maxFractionToDecrease;
		_precisionInConcentration = precisionInConcentration;
		
		// initialize
		setValue(_inputConcentration);
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		
		// reduce concentration
		float newCBulk = getValue() + tstep
			* rate(getCurrentGlobalRateFromMassBalance());
		
		// unless nutrient is below the limit
		if (newCBulk <= _precisionInConcentration) {
			newCBulk = _precisionInConcentration;
		}
		
		setValue(newCBulk);
		System.out.println(newCBulk);
	}

	/**
	 * converts rate into an amount to compute chang in bulk concentration
	 * 
	 * @param r
	 *            the rate in [MT^-1(L of computational volume)^-3]
	 * @return the global rate for this solute species [M(L of reactor)^-3T^-1]
	 */
	private float rate(float r) {
		float a = Model.model().getComputationalVolumeMultiplier()
				/ Model.model().getReactorVolume();
		return a * r;
	}
	
	/**
	 * Guaranty that bulk concentration does not change in steps greater than
	 * FRACTION of the value in previous iteration.
	 * 
	 * @return the maximum time step to fullfil covergence condition
	 */
	private float getMaximumTimeStepForMassBalance() {
		float cBulk = getValue();
		if (cBulk <= _precisionInConcentration)
			// no restrictions are applied if concentration if below the
			// precision
			// threshold
			return Float.POSITIVE_INFINITY;
		float r = rate(getCurrentGlobalRateFromDiffusionReaction());
		// compute mass balance
		// condition 1 - concentration should not change in steps greater than
		// FRACTION
		float t1 = _maxFractionToDecrease * cBulk / Math.abs(r);
		// if t1 is 0 (concentration 0) no time restrain (return infinite)
		return t1;
	}
}