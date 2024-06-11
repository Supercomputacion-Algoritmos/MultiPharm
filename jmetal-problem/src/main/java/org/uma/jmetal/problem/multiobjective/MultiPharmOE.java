package org.uma.jmetal.problem.multiobjective;

import org.uma.jmetal.problem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.List;

// OpenEye imports
import openeye.oechem.*;
import openeye.oezap.*;
import openeye.oeshape.*;


/**
 * Class representing problem MultiPharm
 */
@SuppressWarnings("serial")
public class MultiPharmOE extends AbstractDoubleProblem {

    public static String queryPath;
    public static String targetPath;
    OEGraphMol query = null;
    OEGraphMol target = null;
    double[] coordsTarget = null;
    String queryName = "";
    String targetName = "";
    int contador = 0;
    OEET et = null; 
    OEOverlapPrep prep = null;
    OEExactShapeFunc shapeFunc = null; 
    OEOverlapResults res = null; 


    /**
     * Constructor
     */
    public MultiPharmOE() {

        //String queryPath = "DB01352.mol2";

        //String targetPath = "DB00306.mol2";
        setNumberOfVariables(10);
        setNumberOfObjectives(2);
        setNumberOfConstraints(0);
        setName("MultiPharm");

        // Reading and configuring query molecule
        oemolistream ifs = new oemolistream();
        ifs.SetFormat(OEFormat.MOL2);

        this.query = new OEGraphMol();
        if (!ifs.open(queryPath))
            oechem.OEThrow.Fatal("Unable to open query file.");
        oechem.OEReadMolecule(ifs, query);

        // Align molecule

        // Reading and configuring target molecule
        this.target = new OEGraphMol();
        if (!ifs.open(targetPath))
            oechem.OEThrow.Fatal("Unable to open target file.");
        oechem.OEReadMolecule(ifs, target);
        ifs.close();

        // Shape
        this.prep = new OEOverlapPrep();
        this.prep.SetUseHydrogens(false);
        this.prep.Prep(this.query);
        this.prep.Prep(this.target);

        // ZAP's Stuff
        //query
        oechem.OEAssignBondiVdWRadii(query);
        oechem.OEMMFFAtomTypes(query);
        oechem.OEMMFF94PartialCharges(query);
        // target
        oechem.OEAssignBondiVdWRadii(target);
        oechem.OEMMFFAtomTypes(target);
        oechem.OEMMFF94PartialCharges(target);

        // Get appropriate function to calculate exact shape
        this.shapeFunc = new OEExactShapeFunc();
        this.shapeFunc.SetupRef(this.query);

        
        this.res = new OEOverlapResults();

        // Electrostatic
        this.et = new OEET();
        et.SetRefMol(this.query);

        // Bounds
        List<Double> lowerLimits = new ArrayList<Double>(getNumberOfVariables());
        List<Double> upperLimits = new ArrayList<Double>(getNumberOfVariables());

        this.defineBounds(lowerLimits, upperLimits);

        setLowerLimit(lowerLimits);
        setUpperLimit(upperLimits);

        //float[] xyz = new float[3];

        this.coordsTarget = new double[this.target.NumAtoms() * 3];
        int i3 = 0;
        float[] xyz = new float[3];
        
        for (OEAtomBase atom : this.target.GetAtoms()) {
            this.target.GetCoords(atom, xyz);
            this.coordsTarget[i3] = xyz[0];
            this.coordsTarget[i3 + 1] = xyz[1];
            this.coordsTarget[i3 + 2] = xyz[2];
            i3 += 3;
        }

        this.queryName = this.query.GetTitle();
        this.targetName = this.target.GetTitle();
        this.contador = 0;


        // PRINT THE TWO COMPOUNDS
		/*oemolostream ofs = new oemolostream();
		ofs.open("QUERY.mol2");
		ofs.SetFormat(OEFormat.MOL2);
		oechem.OEWriteMolecule(ofs, this.query);

		ofs = new oemolostream();
		ofs.open("TARGET.mol2");
		ofs.SetFormat(OEFormat.MOL2);
		oechem.OEWriteMolecule(ofs, this.target);
        */

    }

    /**
     * Evaluate() method
     */
    @Override
    public void evaluate(DoubleSolution solution) {
        int numberOfVariables = getNumberOfVariables();
        double[] coords = null;
        double shape = 0;
        double electrostatic = 0;
        
        // CREATE A COPY OF THE TARGET (MIRIAM)
        //OEGraphMol targetcopy = new OEGraphMol(target);
        

        // ROTATE
        //System.out.println(this.contador+ " Parameters: " + this.coordsTarget[0]+ ", "+ this.coordsTarget[1]+ ", " + this.coordsTarget[2] +", "+ this.coordsTarget[3]+", "+ this.coordsTarget[4]+", "+ this.coordsTarget[5]+", "+ this.coordsTarget[6]+", "+ this.coordsTarget[7]+", "+ this.coordsTarget[8]+", "+ this.coordsTarget[9]);
        coords = this.rotate1Axis(this.coordsTarget, solution.getVariableValue(0), solution.getVariableValue(1), solution.getVariableValue(2), solution.getVariableValue(3), solution.getVariableValue(4), solution.getVariableValue(5), solution.getVariableValue(6));

        // DISPLACEMENT
        coords = this.displacement(coords, solution.getVariableValue(7), solution.getVariableValue(8),
                solution.getVariableValue(9));
        // Asign new coords
        this.target.SetCoords(coords);
        // MIRIAM
        //targetcopy.SetCoords(coords);

        /*
		oemolostream ofs = new oemolostream();
		ofs.open("query_"+this.contador+".mol2");
		ofs.SetFormat(OEFormat.MOL2);
		oechem.OEWriteMolecule(ofs, this.query);

		ofs = new oemolostream();
		ofs.open("target_"+this.contador+".mol2");
		ofs.SetFormat(OEFormat.MOL2);
		oechem.OEWriteMolecule(ofs, this.target);
        */

        this.shapeFunc.Overlap(this.target, this.res);
        //MIRIAM
        //this.shapeFunc.Overlap(targetcopy, this.res);
        //System.out.println("overlap");

        //System.out.println(this.contador+ " [" + -this.res.GetTanimoto() + ", " +-this.et.Tanimoto(this.target)+"] "+ solution.getVariableValue(0) + ", "+ solution.getVariableValue(1) + ", "+ solution.getVariableValue(2) + ", "+ solution.getVariableValue(3) + ", "+ solution.getVariableValue(4) + ", "+ solution.getVariableValue(5) + ", "+ solution.getVariableValue(6) + ", "+ solution.getVariableValue(7) + ", "+ solution.getVariableValue(8) + ", "+ solution.getVariableValue(9));
        //this.contador += 1;
        
        
        solution.setObjective(0, -this.res.GetTanimoto()); // obtiene el valor de linea 171
        solution.setObjective(1, -this.et.Tanimoto(this.target));
        //System.out.println("setObj");
        //targetcopy.delete(); 
    }

    private double[] displacement(double[] coords, double deltaX, double deltaY, double deltaZ) {

        int i3 = 0;
        for (int i = 0; i < (coords.length / 3); i++) {
            i3 = i * 3;
            coords[i3] += deltaX;
            coords[i3 + 1] += deltaY;
            coords[i3 + 2] += deltaZ;

        }
        return coords;
    }

    private double[] rotate1Axis(double[] coords, double theta, double x1, double y1, double z1, double x2, double y2,
                                 double z2) {

        double[] new_coords = coords.clone();
        double vectorX = x2 - x1;
        double vectorY = y2 - y1;
        double vectorZ = z2 - z1;

        double angle = theta * 0.5;
        double sinTheta = Math.sin(angle);
        double cosTheta = Math.cos(angle);

        double moduleVector = 0;
        if (vectorX == 0 && vectorY == 0 && vectorZ == 0)
            moduleVector = 1;
        else
            moduleVector = Math.sqrt(vectorX * vectorX + vectorY * vectorY + vectorZ * vectorZ);

        double q1w, q1x, q1y, q1z, q1wConjugated, q1xConjugated, q1yConjugated, q1zConjugated;

        q1w = cosTheta;
        q1x = sinTheta * vectorX / moduleVector;
        q1y = sinTheta * vectorY / moduleVector;
        q1z = sinTheta * vectorZ / moduleVector;

        q1wConjugated = q1w;
        q1xConjugated = -1 * q1x;
        q1yConjugated = -1 * q1y;
        q1zConjugated = -1 * q1z;

        int i3 = 0;

        double atomPositionMinusAW, atomPositionMinusAX, atomPositionMinusAY, atomPositionMinusAZ, part2w, part2x,
                part2y, part2z, part3w, part3x, part3y, part3z;

        for (int i = 0; i < (new_coords.length / 3); i++) {
            i3 = i * 3;
            atomPositionMinusAW = 0;
            atomPositionMinusAX = coords[i3] - x1;
            atomPositionMinusAY = coords[i3 + 1] - y1;
            atomPositionMinusAZ = coords[i3 + 2] - z1;
            part2w = q1w * atomPositionMinusAW - q1x * atomPositionMinusAX - q1y * atomPositionMinusAY
                    - q1z * atomPositionMinusAZ;
            part2x = q1w * atomPositionMinusAX + q1x * atomPositionMinusAW + q1y * atomPositionMinusAZ
                    - q1z * atomPositionMinusAY;
            part2y = q1w * atomPositionMinusAY - q1x * atomPositionMinusAZ + q1y * atomPositionMinusAW
                    + q1z * atomPositionMinusAX;
            part2z = q1w * atomPositionMinusAZ + q1x * atomPositionMinusAY - q1y * atomPositionMinusAX
                    + q1z * atomPositionMinusAW;

            part3w = part2w * q1wConjugated - part2x * q1xConjugated - part2y * q1yConjugated - part2z * q1zConjugated;
            part3x = part2w * q1xConjugated + part2x * q1wConjugated + part2y * q1zConjugated - part2z * q1yConjugated;
            part3y = part2w * q1yConjugated - part2x * q1zConjugated + part2y * q1wConjugated + part2z * q1xConjugated;
            part3z = part2w * q1zConjugated + part2x * q1yConjugated - part2y * q1xConjugated + part2z * q1wConjugated;

            new_coords[i3] = x1 + part3x;
            new_coords[i3 + 1] = y1 + part3y;
            new_coords[i3 + 2] = z1 + part3z;
        }
        return new_coords;
    }

    private void defineBounds(List<Double> lowerLimits, List<Double> upperLimits) {

        lowerLimits.add(0.0);
        upperLimits.add(6.28318530718);

        double x1LowTarget = 0.0;
        double x1UpTarget = 0.0;
        double y1LowTarget = 0.0;
        double y1UpTarget = 0.0;
        double z1LowTarget = 0.0;
        double z1UpTarget = 0.0;

        double xyz[] = new double[3];


        for (OEAtomBase atom : this.target.GetAtoms()) {
            this.target.GetCoords(atom, xyz);
            x1LowTarget = xyz[0];
            x1UpTarget = xyz[0];
            y1LowTarget = xyz[1];
            y1UpTarget = xyz[1];
            z1LowTarget = xyz[2];
            z1UpTarget = xyz[2];
            break;
        }


        for (OEAtomBase atom : this.target.GetAtoms()) {
            this.target.GetCoords(atom, xyz);
            if (xyz[0] < x1LowTarget) {
                x1LowTarget = xyz[0];
            }
            if (xyz[0] > x1UpTarget) {
                x1UpTarget = xyz[0];
            }
            if (xyz[1] < y1LowTarget) {
                y1LowTarget = xyz[1];
            }
            if (xyz[1] > y1UpTarget) {
                y1UpTarget = xyz[1];
            }
            if (xyz[2] < z1LowTarget) {
                z1LowTarget = xyz[2];
            }
            if (xyz[2] > z1UpTarget) {
                z1UpTarget = xyz[2];
            }
        }
        lowerLimits.add(x1LowTarget);
        lowerLimits.add(y1LowTarget);
        lowerLimits.add(z1LowTarget);
        lowerLimits.add(x1LowTarget);
        lowerLimits.add(y1LowTarget);
        lowerLimits.add(z1LowTarget);

        upperLimits.add(x1UpTarget);
        upperLimits.add(y1UpTarget);
        upperLimits.add(z1UpTarget);
        upperLimits.add(x1UpTarget);
        upperLimits.add(y1UpTarget);
        upperLimits.add(z1UpTarget);

        // Delta displacement
        double x1LowQuery = 0.0;
        double x1UpQuery = 0.0;
        double y1LowQuery = 0.0;
        double y1UpQuery = 0.0;
        double z1LowQuery = 0.0;
        double z1UpQuery = 0.0;
        for (OEAtomBase atom : this.query.GetAtoms()) {
            this.query.GetCoords(atom, xyz);
            x1LowQuery = xyz[0];
            x1UpQuery = xyz[0];
            y1LowQuery = xyz[1];
            y1UpQuery = xyz[1];
            z1LowQuery = xyz[2];
            z1UpQuery = xyz[2];
            break;
        }

        for (OEAtomBase atom : this.query.GetAtoms()) {
            this.query.GetCoords(atom, xyz);
            if (xyz[0] < x1LowQuery) {
                x1LowQuery = xyz[0];
            }
            if (xyz[0] > x1UpQuery) {
                x1UpQuery = xyz[0];
            }
            if (xyz[1] < y1LowQuery) {
                y1LowQuery = xyz[1];
            }
            if (xyz[1] > y1UpQuery) {
                y1UpQuery = xyz[1];
            }
            if (xyz[2] < z1LowQuery) {
                z1LowQuery = xyz[2];
            }
            if (xyz[2] > z1UpQuery) {
                z1UpQuery = xyz[2];
            }
        }

        double minDeltaX = x1LowTarget - x1LowQuery;
        double minDeltaY = y1LowTarget - y1LowQuery;
        double minDeltaZ = z1LowTarget - z1LowQuery;
        double maxDeltaX = x1UpTarget - x1UpQuery;
        double maxDeltaY = y1UpTarget - y1UpQuery;
        double maxDeltaZ = z1UpTarget - z1UpQuery;

        double deltaX = Math.max(Math.abs(minDeltaX), Math.abs(maxDeltaX));
        double deltaY = Math.max(Math.abs(minDeltaY), Math.abs(maxDeltaY));
        double deltaZ = Math.max(Math.abs(minDeltaZ), Math.abs(maxDeltaZ));

        lowerLimits.add(-Math.abs(deltaX));
        lowerLimits.add(-Math.abs(deltaY));
        lowerLimits.add(-Math.abs(deltaZ));

        upperLimits.add(Math.abs(deltaX));
        upperLimits.add(Math.abs(deltaY));
        upperLimits.add(Math.abs(deltaZ));
    }

    private void centerMol(OEGraphMol mol) {
        double centroidX = 0.0, centroidY = 0.0, centroidZ = 0.0;
        int count = 0;

        float xyz[] = new float[3];
        for (OEAtomBase atom : mol.GetAtoms()) {
            this.query.GetCoords(atom, xyz);
            centroidX += xyz[0];
            centroidY += xyz[1];
            centroidZ += xyz[2];
            count += 1;
        }
        centroidX /= count;
        centroidY /= count;
        centroidZ /= count;

        double[] newCoords = new double[3];
        for (OEAtomBase atom : mol.GetAtoms()) {
            mol.GetCoords(atom, xyz);
            newCoords[0] = xyz[0] - centroidX;
            newCoords[1] = xyz[1] - centroidY;
            newCoords[2] = xyz[2] - centroidZ;
            mol.SetCoords(atom, newCoords);
        }
    }

    private void pcaProjectMol(OEGraphMol mol) {

    }

}
