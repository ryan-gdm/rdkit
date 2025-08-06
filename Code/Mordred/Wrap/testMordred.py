from rdkit import Chem
import numpy as np
import rdMordred as rd
import time
import pandas as pd


def CalcMordred(mol: Chem.Mol) -> np.ndarray | None:
    version = 2
    doExEstate = True


    descriptor_funcs = [
        (rd.CalcABCIndex, [mol]),
        (rd.CalcAcidBase, [mol]),
        (rd.CalcAdjacencyMatrix, [mol, version]),
        (rd.CalcAromatic, [mol]),
        (rd.CalcAtomCount, [mol, version]),
        (rd.CalcAutocorrelation, [mol]),
        (rd.CalcBCUT, [mol]),
        (rd.CalcBalabanJ, [mol]),
        (rd.CalcBaryszMatrix, [mol]),
        (rd.CalcBertzCT, [mol]),
        (rd.CalcBondCount, [mol]),
        (rd.CalcRNCGRPCG, [mol]),
        (rd.CalcCarbonTypes, [mol, version]),
        (rd.CalcChi, [mol]),
        (rd.CalcConstitutional, [mol]),
        (rd.CalcDetourMatrix, [mol]),
        (rd.CalcDistanceMatrix, [mol, version]),
        (rd.CalcEState, [mol, doExEstate]),
        (rd.CalcEccentricConnectivityIndex, [mol]),
        (rd.CalcExtendedTopochemicalAtom, [mol]),
        (rd.CalcFragmentComplexity, [mol]),
        (rd.CalcFramework, [mol]),
        (rd.CalcHydrogenBond, [mol]),
        (rd.CalcLogS, [mol]),
        (rd.CalcInformationContent, [mol, 5]),
        (rd.CalcKappaShapeIndex, [mol]),
        (rd.CalcLipinski, [mol]),
        (rd.CalcMcGowanVolume, [mol]),
        (rd.CalcMoeType, [mol]),
        (rd.CalcMolecularDistanceEdge, [mol]),
        (rd.CalcMolecularId, [mol]),
        (rd.CalcPathCount, [mol]),
        (rd.CalcPolarizability, [mol]),
        (rd.CalcRingCount, [mol]),
        (rd.CalcRotatableBond, [mol]),
        (rd.CalcSLogP, [mol]),
        (rd.CalcTopoPSA, [mol]),
        (rd.CalcTopologicalCharge, [mol]),
        (rd.CalcTopologicalIndex, [mol]),
        (rd.CalcVdwVolumeABC, [mol]),
        (rd.CalcVertexAdjacencyInformation, [mol]),
        (rd.CalcWalkCount, [mol]),
        (rd.CalcWeight, [mol]),
        (rd.CalcWienerIndex, [mol]),
        (rd.CalcZagrebIndex, [mol]),
        (rd.CalcPol, [mol]),
        (rd.CalcMR, [mol]),
        (rd.CalcODT, [mol]),
        (rd.CalcFlexibility, [mol]),
        (rd.CalcSchultz, [mol]),
        (rd.CalcAlphaKappaShapeIndex, [mol]),
        (rd.CalcHEState, [mol]),
        (rd.CalcBEState, [mol]),
        (rd.CalcAbrahams, [mol]),
        (rd.CalcANMat, [mol]),
        (rd.CalcASMat, [mol]),
        (rd.CalcAZMat, [mol]),
        (rd.CalcDSMat, [mol]),
        (rd.CalcDN2Mat, [mol]),
        (rd.CalcFrags, [mol]),
        (rd.CalcAddFeatures, [mol])
    ]

    try:
        results = [np.atleast_1d(func(*params)) for func, params in descriptor_funcs]
        return np.concatenate(results)
    except Exception as e:
        print(f"Error processing molecule {mol}: {e}")
        return None




def benchmark_calc_function(mols):
    print("\nBenchmarking: CalcMordred(smiles) per molecule")
    times = []
    results = []

    for mol in mols:
        start = time.time()
        result = CalcMordred(mol)
        end = time.time()

        times.append(end - start)
        results.append(result)

    times = np.array(times)
    print(f"Processed {len(results)} molecules")
    print(f"Total time: {times.sum():.2f} seconds")
    print(f"Avg time per mol: {times.mean():.6f} s | Min: {times.min():.6f} s | Max: {times.max():.6f} s | Std: {times.std():.6f} s")
    print(f"Example output shape: {results[0].shape if hasattr(results[0], 'shape') else type(results[0])}")

    return times

def benchmark_calculator_class(mols):
    print("\nBenchmarking: Calculator.map(mol) per molecule")
    calc = Calculator(descriptors)
    times = []
    results = []

    for mol in mols:
        start = time.time()
        result = list(calc.map([mol]))[0]
        end = time.time()

        times.append(end - start)
        results.append(result)

    times = np.array(times)
    print(f"Processed {len(results)} molecules")
    print(f"Total time: {times.sum():.2f} seconds")
    print(f"Avg time per mol: {times.mean():.6f} s | Min: {times.min():.6f} s | Max: {times.max():.6f} s | Std: {times.std():.6f} s")

    descriptor_names = [str(d) for d in calc.descriptors]
    df = pd.DataFrame(results, columns=descriptor_names)
    print(f"Number of descriptors: {df.shape[1]}")

    return times

if __name__ == "__main__":
    df = pd.read_csv("data/example.csv")
    smiles_list = df["SMILES"].dropna().tolist()
    mols = [Chem.MolFromSmiles(s) for s in smiles_list if Chem.MolFromSmiles(s)]

    print(f"Total molecules to process: {len(mols)}")

    times_custom = benchmark_calc_function(mols)
    times_mordred = benchmark_calculator_class(mols)

    print("\nSummary:")
    print(f"Custom CalcMordred avg time: {np.mean(times_custom)*1000:.2f} ms/molecule")
    print(f"Mordred Calculator avg time  : {np.mean(times_mordred)*1000:.2f} ms/molecule")
