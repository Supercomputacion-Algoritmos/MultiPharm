package org.uma.jmetal.problem.multiobjective;

import openeye.oechem.*;
import openeye.oeshape.OEExactShapeFunc;
import openeye.oeshape.OEOverlapPrep;
import openeye.oeshape.OEOverlapResults;
import openeye.oezap.OEET;
import org.uma.jmetal.problem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// OpenEye imports

/**
 * Class representing problem MultiPharm
 */
@SuppressWarnings("serial")
public class MultiPharm extends AbstractDoubleProblem {

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
    boolean sameVanDerWaalsRadius = true;
    double [] weightQuery = null;
    double [] radiusQuery = null;
    double [] weightTarget= null;
    double [] radiusTarget = null;
    double [] coordsQuery = null;
    double tanimotoQuery = 0;
    double tanimotoTarget = 0;

    Map<String, Double> defaultVDWR = new HashMap<String, Double>() {{
        put("", 0.0);
        put("H", 1.20);
        put("He", 1.40);
        put("Li", 1.82);
        put("Be", 2.00);
        put("B", 2.00);
        put("C", 1.70);
        put("N", 1.55);
        put("O", 1.52);
        put("F", 1.47);
        put("Ne", 1.54);
        put("Na", 2.27);
        put("Mg", 1.73);
        put("Al", 2.00);
        put("Si", 2.10);
        put("P", 1.80);
        put("S", 1.80);
        put("Cl", 1.75);
        put("Ar", 1.88);
        put("K", 2.75);
        put("Ca", 2.00 );
        put("Sc", 2.00);
        put("Ti", 2.00);
        put("V", 2.00);
        put("Cr", 2.00);
        put("Mn", 2.00);
        put("Fe", 2.00);
        put("Co", 2.00);
        put("Ni", 1.63);
        put("Cu", 1.40);
        put("Zn", 1.39);
        put("Ga", 1.87);
        put("Ge", 2.00);
        put("As", 1.85);
        put("Se", 1.90);
        put("Br", 1.85 );
        put("Kr", 2.02);
        put("Rb", 2.00);
        put("Sr", 2.00);
        put("Y", 2.00);
        put("Zr", 2.00);
        put("Nb", 2.00);
        put("Mo", 2.00);
        put("Tc", 2.00);
        put("Ru", 2.00 );
        put("Rh", 2.00);
        put("Pd", 1.63);
        put("Ag", 1.72);
        put("Cd", 1.58);
        put("In", 1.93);
        put("Sn", 2.17);
        put("Sb", 2.00);
        put("Te", 2.06);
        put("I", 1.98);
        put("Xe", 2.16 );
        put("Cs", 2.00);
        put("Ba", 2.00 );
        put("La", 2.00);
        put("Ce", 2.00);
        put("Pr", 2.00);
        put("Nd", 2.00);
        put("Pm", 2.00);
        put("Sm", 2.00);
        put("Eu", 2.00);
        put("Gd", 2.00);
        put("Tb", 2.00);
        put("Dy", 2.00 );
        put("Ho", 2.00);
        put("Er", 2.00 );
        put("Tm", 2.00);
        put("Yb", 2.00);
        put("Lu", 2.00);
        put("Hf", 2.00 );
        put("Ta", 2.00);
        put("W", 2.00);
        put("Re", 2.00);
        put("Os", 2.00);
        put("Ir", 2.00 );
        put("Pt", 1.72);
        put("Au", 1.66);
        put("Hg", 1.55);
        put("Tl", 1.96);
        put("Pb", 2.02);
        put("Bi", 2.00);
        put("Po", 2.00);
        put("At", 2.00);
        put("Rn", 2.00);
        put("Fr", 2.00);
        put("Ra", 2.00);
        put("Ac", 2.00);
        put("Th", 2.00);
        put("Pa", 2.00);
        put("U", 1.86);
        put("Np", 2.00);
        put("Pu", 2.00);
        put("Am", 2.00);
        put("Cm", 2.00);
        put("Bk", 2.00);
        put("Cf", 2.00);
        put("Es", 2.00);
        put("Fm", 2.00);
        put("Md", 2.00);
        put("No", 2.00);
        put("Lr", 2.00);
        put("Rf", 2.00 );
        put("Db", 2.00 );
        put("Sg", 2.00);
        put("Bh", 2.00);
        put("Hs", 2.00);
        put("Mt", 2.00);
        put("Ds", 2.00);
        put("Rg", 2.00);
    }};

    /**
     * Constructor
     */
    public MultiPharm() {

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

        // ZAP's Stuff
        //query
        oechem.OEAssignBondiVdWRadii(query);
        oechem.OEMMFFAtomTypes(query);
        oechem.OEMMFF94PartialCharges(query);
        // target
        oechem.OEAssignBondiVdWRadii(target);
        oechem.OEMMFFAtomTypes(target);
        oechem.OEMMFF94PartialCharges(target);

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

        // WEIGHT and RADIUS QUERY
        this.weightQuery = new double[this.query.NumAtoms()];
        this.radiusQuery = new double[this.query.NumAtoms()];
        int counter = 0;
        for (OEAtomBase atom : this.query.GetAtoms()) {
            this.radiusQuery[counter] = 1.80;
            counter +=1 ;
        }

        counter = 0;
        for (OEAtomBase atom : this.query.GetAtoms()) {
            this.query.GetCoords(atom, xyz);
            this.weightQuery[counter] = calculateWeightWEGA(xyz[0],xyz[1],xyz[2],this.query, this.radiusQuery,counter,sameVanDerWaalsRadius);
            counter +=1;

        }

        // WEIGHT and RADIUS TARGET
        this.weightTarget = new double[this.target.NumAtoms()];
        this.radiusTarget = new double[this.target.NumAtoms()];
        counter = 0;
        for (OEAtomBase atom : this.target.GetAtoms()) {
            this.radiusTarget[counter] = 1.8;//defaultVDWR.get(atom.GetName());
            counter += 1;
        }


        counter = 0;
        for (OEAtomBase atom : this.target.GetAtoms()) {
            this.target.GetCoords(atom, xyz);
            this.weightTarget[counter] = calculateWeightWEGA(xyz[0],xyz[1],xyz[2],this.target, this.radiusTarget,counter,sameVanDerWaalsRadius);
            counter +=1;

        }

        this.queryName = this.query.GetTitle();
        this.targetName = this.target.GetTitle();
        this.contador = 0;
        this.coordsQuery = new double[this.query.NumAtoms()*3];
        this.query.GetCoords(this.coordsQuery);

        this.tanimotoQuery = overlapWEGA(this.coordsQuery,this.weightQuery,this.radiusQuery,this.coordsQuery,this.weightQuery,this.radiusQuery,this.sameVanDerWaalsRadius);
        double [] tempCoordsTarget = new double[this.target.NumAtoms()*3];
        this.target.GetCoords(tempCoordsTarget);
        this.tanimotoTarget = overlapWEGA(tempCoordsTarget,this.weightTarget,this.radiusTarget,tempCoordsTarget,this.weightTarget,this.radiusTarget, this.sameVanDerWaalsRadius);

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
        
        if ((Math.abs(solution.getVariableValue(1)-solution.getVariableValue(4))<1.0E-16)&&(Math.abs(solution.getVariableValue(2)-solution.getVariableValue(5))<1.0E-16)&&(Math.abs(solution.getVariableValue(3)-solution.getVariableValue(6))<1.0E-16)){
            //System.out.println("IF");
            solution.setObjective(0, 0.0);
            solution.setObjective(1, 0.0);
        }
        else{
            // ROTATE
            //System.out.println(this.contador+ " Parameters: " + this.coordsTarget[0]+ ", "+ this.coordsTarget[1]+ ", " + this.coordsTarget[2] +", "+ this.coordsTarget[3]+", "+ this.coordsTarget[4]+", "+ this.coordsTarget[5]+", "+ this.coordsTarget[6]+", "+ this.coordsTarget[7]+", "+ this.coordsTarget[8]+", "+ this.coordsTarget[9]);
            coords = this.rotate1Axis(this.coordsTarget, solution.getVariableValue(0), solution.getVariableValue(1), solution.getVariableValue(2), solution.getVariableValue(3), solution.getVariableValue(4), solution.getVariableValue(5), solution.getVariableValue(6));
    
    
            // DISPLACEMENT
            coords = this.displacement(coords, solution.getVariableValue(7), solution.getVariableValue(8),
                    solution.getVariableValue(9));
            // Asign new coords
            this.target.SetCoords(coords);
    
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
    
            double VAB = overlapWEGA(this.coordsQuery,this.weightQuery,this.radiusQuery,coords,this.weightTarget,this.radiusTarget,this.sameVanDerWaalsRadius);
            double shapeTanimoto = (VAB/(this.tanimotoQuery+this.tanimotoTarget-VAB));
            //System.out.println(this.contador+ " [" + -shapeTanimoto + ", " +-this.et.Tanimoto(this.target)+"] "+ solution.getVariableValue(0) + ", "+ solution.getVariableValue(1) + ", "+ solution.getVariableValue(2) + ", "+ solution.getVariableValue(3) + ", "+ solution.getVariableValue(4) + ", "+ solution.getVariableValue(5) + ", "+ solution.getVariableValue(6) + ", "+ solution.getVariableValue(7) + ", "+ solution.getVariableValue(8) + ", "+ solution.getVariableValue(9));
    
            //System.out.println(this.contador+ " [" + -this.res.GetTanimoto() + ", " +-this.et.Tanimoto(this.target)+"] "+ solution.getVariableValue(0) + ", "+ solution.getVariableValue(1) + ", "+ solution.getVariableValue(2) + ", "+ solution.getVariableValue(3) + ", "+ solution.getVariableValue(4) + ", "+ solution.getVariableValue(5) + ", "+ solution.getVariableValue(6) + ", "+ solution.getVariableValue(7) + ", "+ solution.getVariableValue(8) + ", "+ solution.getVariableValue(9));
            //this.contador += 1;
            solution.setObjective(0, -shapeTanimoto);
            solution.setObjective(1, -this.et.Tanimoto(this.target));
        }
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
        
        x1LowTarget = Math.min(x1LowTarget, -0.5);
        y1LowTarget = Math.min(y1LowTarget, -0.5);
        z1LowTarget = Math.min(z1LowTarget, -0.5);
        
        x1UpTarget = Math.max(x1UpTarget, 0.5);
        y1UpTarget = Math.max(y1UpTarget, 0.5);
        z1UpTarget = Math.max(z1UpTarget, 0.5);
        
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
        
        deltaX = Math.max(deltaX, 0.5);
        deltaY = Math.max(deltaY, 0.5);
        deltaZ = Math.max(deltaZ, 0.5);

        lowerLimits.add(-Math.abs(deltaX));
        lowerLimits.add(-Math.abs(deltaY));
        lowerLimits.add(-Math.abs(deltaZ));

        upperLimits.add(Math.abs(deltaX));
        upperLimits.add(Math.abs(deltaY));
        upperLimits.add(Math.abs(deltaZ));
        
        /*System.out.println("lowerLimits");
        for (int i =0; i<10; i++){
            System.out.println(lowerLimits.get(i));
        }
        System.out.println("upperLimits");
        for (int i =0; i<10; i++){
            System.out.println(upperLimits.get(i));
            }*/
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

    private double overlapWEGA(double[] queryAtomsXYZ, double[] queryWeightAtoms, double[] queryRadiusAtoms,double[] targetAtomsXYZ, double[] targetWeightAtoms, double[] targetRadiusAtoms,boolean sameVanDerWaalsRadius){
        double overlap = 0.0;
        for (int i = 0; i< (queryAtomsXYZ.length)/3;i++){
            int iteradori = i*3;
            for (int j=0;j< (targetAtomsXYZ.length)/3;j++){
                int iteradorj = j * 3;
                double Vij = calculateOverlapVolumeAtomsWEGA(queryAtomsXYZ[iteradori],
                        queryAtomsXYZ[iteradori + 1], queryAtomsXYZ[iteradori + 2],
                        targetAtomsXYZ[iteradorj], targetAtomsXYZ[iteradorj + 1],
                        targetAtomsXYZ[iteradorj + 2], sameVanDerWaalsRadius,
                        queryRadiusAtoms[i], targetRadiusAtoms[j]);

                double WeightVij = queryWeightAtoms[i] * targetWeightAtoms[j] * Vij;
                overlap += WeightVij;
            }
        }
        return overlap;
    }

    private double calculateOverlapVolumeAtomsWEGA(double atom1X,
                                                   double atom1Y, double atom1Z, double atom2X, double atom2Y,
                                                   double atom2Z, boolean sameVanDerWaalsRadius, double radiusFixed,
                                                   double radiusVariable) {

        double Vij = 0;
        double tempX = atom1X - atom2X;
        double tempY = atom1Y - atom2Y;
        double tempZ = atom1Z - atom2Z;
        double Rij2 = tempX * tempX + tempY * tempY + tempZ * tempZ;

        if (sameVanDerWaalsRadius) {
            //		double Kij2 = exp(Rij2 * -0.3731438999881213);
            double Kij = Math.exp(Rij2 * -0.3731438999881213);

            //Vij = 7.99984656 * Kij * pow(2.104813085306902, 1.5);
            Vij = 24.428790199 * Kij;
        } else {
            double pi = 3.14159265358;
            double alphai = 2.417972471923026 / (radiusFixed * radiusFixed);
            double alphaj = 2.417972471923026 / (radiusVariable * radiusVariable);

            double Kij = Math.exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
            Vij = 7.99984656 * Kij * Math.pow((pi / (alphai + alphaj)), 1.5);

        }

        return Vij;

    }

    private double calculateWeightWEGA(double atomX, double atomY,
                                       double atomZ, OEGraphMol molecule, double []radius, int iAtom,
                                       boolean sameVanDerWaalsRadius){
        double pi = 3.14159265358;
        double vi = 4 * pi * Math.pow(radius[iAtom], 3) / 3;

        //double p = 2.8284;

        double vij = 0;
        double k = 0.8665;
        int iterador = 0;
        double [] xyz = new double[3];
        for (OEAtomBase atom : molecule.GetAtoms()) {

            if (iterador == iAtom) {
                iterador+=1;
                continue;
            }
            molecule.GetCoords(atom, xyz);
            vij += calculateOverlapVolumeAtomsWEGA(atomX, atomY, atomZ,
                    xyz[0], xyz[1], xyz[2],
                    sameVanDerWaalsRadius, radius[iAtom], radius[iterador]);
            //System.out.println("VIJ de "+iAtom+ ": " + vij);
            iterador+=1;
        }

        double result = vi / (vi + k * vij);

        return result;
    }
    private void pcaProjectMol(OEGraphMol mol) {

    }

}
