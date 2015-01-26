package nl.tudelft.bt.model;

import java.io.Serializable;

/**
 * Implements 3D vector of continuous spatial coordinates
 * 
 * @author João Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ContinuousCoordinate implements Serializable {
	public float x;
	public float y;
	public float z;

	/**
	 * Constructor (creates 0 length vector)
	 */
	public ContinuousCoordinate() {
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	}

	/**
	 * Constructor.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 */
	public ContinuousCoordinate(float x, float y, float z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	/**
	 * set x = y = z = 0
	 */
	public void reset() {
		x = 0.0f;
		y = 0.0f;
		z = 0.0f;
	}
	/**
	 * Print coordinates to string
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return x + ",\t" + y + ",\t" + z;
	}
}