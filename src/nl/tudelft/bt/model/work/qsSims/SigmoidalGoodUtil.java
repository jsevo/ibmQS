/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.qsSims;

import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SigmoidalGoodUtil extends ProcessFactor {
	protected MultigridVariable _publicGood;

	private float _steepness;

	private float _goodHalfPoint;

	/**
	 * @param c
	 * @param k
	 */
	public SigmoidalGoodUtil(MultigridVariable publicGood, float steepness, float goodHalfPoint) {
		_publicGood = publicGood;
		_steepness = steepness;
		_goodHalfPoint = goodHalfPoint;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concPublicGood = _publicGood.getValue();
		return (float) Math.tanh((concPublicGood - _goodHalfPoint) * _steepness) + 1f;
	}

	public float getMaximumValue() {
		float concPublicGood = _publicGood.getMaximumValue();
		return (float) (Math.tanh((concPublicGood - _goodHalfPoint) * _steepness)) + 1f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		float concPublicGood = _publicGood.getValue();
		return (float) (_steepness * (1 - Math.tanh((concPublicGood - _goodHalfPoint) * _steepness)));
	}
}
