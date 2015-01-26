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
public class SoluteLevelCheck extends ProcessFactor {
	protected MultigridVariable _solute;

	/**
	 * @param solute
	 * @param threshold
	 */
	public SoluteLevelCheck(MultigridVariable solute){
		_solute = solute;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concSolute = _solute.getValue();
		if(concSolute != 7f) {for(int i = 0; i < 100; i++)
		System.out.println(concSolute);}
		else System.out.println(concSolute);
			return 1;
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
