/* 
 * Created on 16-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.exceptions;

/**
 * Exception thrown when biofilm is not growing for lack of nutrients
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiofilmProblemException extends ModelRuntimeException {

	/**
	 * 
	 */
	public BiofilmProblemException() {
		super("Biofilm is not growing for lack of nutrients");
	}
	/**
	 * @param message
	 */
	public BiofilmProblemException(String message) {
		super(message);
	}
}
