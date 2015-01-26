/*
 * Created on April 28, 2009
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class StepOrNoStep extends ProcessFactor {
	protected MultigridVariable _solute;
	private float _threshold;
	private boolean _ifStep;

	/**
	 * @param solute
	 * @param threshold
	 */
	public StepOrNoStep(MultigridVariable solute, float threshold, boolean ifStep){
		_solute = solute;
		_threshold = threshold;
		_ifStep = ifStep;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		if (_ifStep)
		{
			float concPublicGood = _solute.getValue();
			if (concPublicGood < 0) System.out.println(concPublicGood);
			if (concPublicGood > _threshold){
				return 1f;
			} else {
				return 0f;			
			}
		}
		else
			return 1f;
	}

	public float getMaximumValue() {
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		return 0;
	}
}
