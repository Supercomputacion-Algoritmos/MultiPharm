from openeye import oeshape
from openeye import oezap
from openeye import oechem
import sys
import math
# if it is not possible to install the library
#sys.path.append("/Users/savins/repositorios/jMetalPy/libs/OpenEye-toolkits-python3-osx-x64-2019.10.2")


def displace(molecule_atoms, deltaX, deltaY, deltaZ):
    """Displace atoms coordinates according to the parameters.

    Args:
        molecule_atoms ([type]): array of the coordinates of the atoms. Each coordinate is in one position so the coordinates of the first atom are the 0, 1 and 2 position of the array.
        deltaX ([type]): displacement in X
        deltaY ([type]): displacement in Y
        deltaZ ([type]): displacement in Z
    """

    counter = 0
    for i in range(0, int(len(molecule_atoms)/3)):
        molecule_atoms[counter] += deltaX
        molecule_atoms[counter+1] += deltaY
        molecule_atoms[counter+2] += deltaZ
        counter += 3

    return molecule_atoms


def rotateAccording1Axis(molecule_atoms, theta, x1, y1, z1, x2, y2, z2):
    vectorx = x2-x1
    vectory = y2-y1
    vectorz = z2-z1

    sinTheta = math.sin(theta/2)
    cosTheta = math.cos(theta/2)

    moduleVector = math.sqrt(vectorx * vectorx + vectory * vectory + vectorz * vectorz)

    q1W = cosTheta
    q1X = sinTheta * vectorx / moduleVector
    q1Y = sinTheta * vectory / moduleVector
    q1Z = sinTheta * vectorz / moduleVector

    q1ConjugatedW = q1W
    q1ConjugatedX = -q1X
    q1ConjugatedY = -q1Y
    q1ConjugatedZ = -q1Z

    counter = 0
    for i in range(0, int(len(molecule_atoms)/3)):
        atomPositionMinuxAX = molecule_atoms[counter] - x1
        atomPositionMinuxAY = molecule_atoms[counter+1] - y1
        atomPositionMinuxAZ = molecule_atoms[counter+2] - z1

        part2W, part2X, part2Y, part2Z = multiply(
            q1W, q1X, q1Y, q1Z, 0, atomPositionMinuxAX, atomPositionMinuxAY, atomPositionMinuxAZ)
        
        part3W, part3X, part3Y, part3Z = multiply(
            part2W, part2X, part2Y, part2Z, q1ConjugatedW, q1ConjugatedX, q1ConjugatedY, q1ConjugatedZ)

        pW, pX, pY, pZ = add(0, x1, y1, z1, part3W, part3X, part3Y, part3Z)

        molecule_atoms[counter] = pX
        molecule_atoms[counter+1] = pY
        molecule_atoms[counter+2] = pZ
        counter += 3
    return molecule_atoms


def add(aW, aX, aY, aZ, bW, bX, bY, bZ):
    qW = aW + bW
    qX = aX + bX
    qY = aY + bY
    qZ = aZ + bZ
    return qW, qX, qY, qZ


def multiply(aW, aX, aY, aZ, bW, bX, bY, bZ):
    qW = aW * bW - aX * bX - aY * bY - aZ * bZ
    qX = aW * bX + aX * bW + aY * bZ - aZ * bY
    qY = aW * bY - aX * bZ + aY * bW + aZ * bX
    qZ = aW * bZ + aX * bY - aY * bX + aZ * bW
    return qW, qX, qY, qZ


if __name__ == "__main__":
    if (len(sys.argv) != 12):
        print("HELP:\nThis script rotates and displaces a molecule which is in mol2 format. The output file has the same name than the input one adding the termination \"_Rotated\".\n\nHow to use it: \n\t python {} <mol2File> <angle> <x1> <y1> <z1> <x2> <y2> <z2> <deltaX> <deltaY> <deltaZ>".format(
            sys.argv[0]))
        sys.exit()
    coords = []

    # Read the molecule
    ifs = oechem.oemolistream()
    ifs.SetFormat(oechem.OEFormat_MOL2)
    #start = time.time()
    # Reading molecules refmol (query)
    query = oechem.OEGraphMol()
    if not ifs.open(sys.argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sys.argv[1])
    oechem.OEReadMolecule(ifs, query)
    for atom in query.GetAtoms():
        x, y, z = query.GetCoords(atom)
        coords.append(x)
        coords.append(y)
        coords.append(z)

    coords2 = rotateAccording1Axis(coords, float(sys.argv[2]), float(
        sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]))

    coords3 = displace(coords2, float(sys.argv[9]), float(
        sys.argv[10]), float(sys.argv[11]))

    query.SetCoords(coords3)
    ofs = oechem.oemolostream(
        "{}_Rotated.mol2".format(sys.argv[1].split(".")[0]))
    ofs.SetFormat(oechem.OEFormat_MOL2)
    oechem.OEWriteMolecule(ofs, query)
