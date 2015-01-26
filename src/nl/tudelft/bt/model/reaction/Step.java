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
public class Step extends ProcessFactor {
	protected MultigridVariable _solute;
	private float _threshold;

	/**
	 * @param solute
	 * @param threshold
	 */
	public Step(MultigridVariable solute, float threshold){
		_solute = solute;
		_threshold = threshold;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concPublicGood = _solute.getValue();
		
		//concPublicGood = 10f;
		if (concPublicGood > _threshold){
			return 1f;
		} else {
			return 0f;			
		}
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
