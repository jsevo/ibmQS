/* 
 * Created on 10-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset.functions;

import nl.tudelft.bt.model.ContinuousCoordinate;

/**
 * Implements constant volumetric detachment
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class VolumetricDetachment
	extends DetachmentSpeedFunction {
	float _detachmentRate;
	
	/**
	 * @param rate
	 */
	public VolumetricDetachment(float rate) {
		_detachmentRate = rate;
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate, float)
	 */
	public float getValue(ContinuousCoordinate c) {
		return _detachmentRate;
	}
	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.detachment.DetachmentSpeedFunction#detachmentIsOff()
	 */
	public boolean detachmentIsOff() {
		return _detachmentRate == 0;
	}
}
