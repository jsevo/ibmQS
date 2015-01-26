package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Everything is inside the boundary layer. Solves the system everywhere.
 * Substrate must come from influx
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class NoBoundaryLayer extends BoundaryLayer {

	/**
	 * @throws MultigridSystemNotSetException
	 */
	public NoBoundaryLayer() throws MultigridSystemNotSetException {
		super();
	}

	public void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		int n = bl.length;
		int m = bl[0].length;
		int l = bl[0][0].length;
	}
}