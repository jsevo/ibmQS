/*
 * File created originally on Jul 7, 2005
 */
package nl.tudelft.bt.model.multigrid.boundary_conditions;

import java.io.Serializable;

/**
 * Implements boundary conditions for a granule type of morphology where bottom
 * all borders are constant value
 * 
 * @author jxavier
 */
public class GranuleBoundaryConditions implements BoundaryConditions,
		Serializable {

	// // ZERO-FLUX BOUNDARIES EVERYWHERE
	// public void refreshBoundaryConditions(float[][][] u) {
	// int l = u[0][0].length - 2;
	// int m = u[0].length - 2;
	// int n = u.length - 2;
	//
	// //zero flux borders everywhere
	// for (int i = 1; i <= n; i++) {
	// for (int j = 1; j <= m; j++) {
	// // sides
	// u[i][j][0] = u[i][j][1];
	// u[i][j][l + 1] = u[i][j][l];
	// }
	// for (int k = 1; k <= l; k++) {
	// // sides
	// u[i][0][k] = u[i][1][k];
	// u[i][m + 1][k] = u[i][m][k];
	// }
	// }
	// for (int j = 1; j <= m; j++) {
	// for (int k = 1; k <= l; k++) {
	// // bottom and top
	// u[0][j][k] = u[1][j][k];
	// u[n+1][j][k] = u[n][j][k];
	// }
	// }
	// }

	// // CYCLIC BOUNDARIES EVERYWHERE
	public void refreshBoundaryConditions(float[][][] u) {
		int l = u[0][0].length - 2;
		int m = u[0].length - 2;
		int n = u.length - 2;

		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				// cyclic borders (sides)
				u[i][j][0] = u[i][j][l];
				u[i][j][l + 1] = u[i][j][1];
			}
			for (int k = 1; k <= l; k++) {
				// cyclic borders (sides)
				u[i][0][k] = u[i][m][k];
				u[i][m + 1][k] = u[i][1][k];
			}
		}
		for (int j = 1; j <= m; j++) {
			for (int k = 1; k <= l; k++) {
				u[0][j][k] = u[n + 1][j][k];
				u[n + 1][j][k] = u[0][j][k];
			}
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#isCarrier(int,
	 *      int, int)
	 */
	public boolean isCarrier(int i, int j, int k) {
		// returns false since the carrier for planar geometry is outrside
		// the computational domain
		return false;
	}
}