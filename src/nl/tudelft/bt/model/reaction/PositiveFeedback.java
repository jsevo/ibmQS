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
public class PositiveFeedback extends ProcessFactor {
	protected MultigridVariable _solute;
	private float _threshold;
	private float _enhancementFactor;

	/**
	 * @param solute
	 * @param threshold
	 */
	public PositiveFeedback(MultigridVariable solute, float threshold, float enhancementFactor){
		_solute = solute;
		_threshold = threshold;
		_enhancementFactor = enhancementFactor;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concSolute = _solute.getValue();
		if (concSolute > _threshold){
			return _enhancementFactor;
		} else {
			return 1;			
		}
	}

	public float getMaximumValue() {
		return _enhancementFactor;
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
