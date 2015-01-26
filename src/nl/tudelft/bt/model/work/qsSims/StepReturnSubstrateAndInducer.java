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
 * @para 
 */
public class StepReturnSubstrateAndInducer extends ProcessFactor {
	protected MultigridVariable _substrate;

	protected MultigridVariable _inducer;

	private float _Ks;

	private float _inducerThreshold;

	/**
	 * @param c
	 * @param k
	 */
	public StepReturnSubstrateAndInducer(MultigridVariable substrate,
			MultigridVariable inducer, float Ks, float inducerThreshold) {
		_substrate = substrate;
		_inducer = inducer;
		_Ks = Ks;
		_inducerThreshold = inducerThreshold;		
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concSubstrate = _substrate.getValue();
		float concInducer = _inducer.getValue();
		//System.out.println("InducerConc:" + _inducer.getValue() + _inducerThreshold);
		//concInducer = 10f;
		if (concInducer < _inducerThreshold) {
		//if(false) {
			return 0;
		} else {
			System.out.println("InducerConc2:" + _inducer.getValue() + _inducerThreshold);
			return concSubstrate / (_Ks + concSubstrate);
		}
	}

	public float getMaximumValue() {
		float concSubstrate = _substrate.getMaximumValue();
		float concInducer = _inducer.getMaximumValue();
		concInducer = (concInducer < 0 ? 0 : concInducer);
		return concSubstrate / (_Ks + concSubstrate);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		float concInducer = _inducer.getValue();
		if (concInducer < _inducerThreshold)
		//if(false)
			return 0;
		if (c == _substrate) {
			float conc1 = _substrate.getValue();
			return _Ks / (_Ks + conc1) / (_Ks + conc1);
		} else if (c == _inducer) {
			return 0;
		}
		return 0f; // derivative for other chemicals
	}
}
