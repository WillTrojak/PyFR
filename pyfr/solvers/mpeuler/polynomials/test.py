from pyfr.solvers.mpeuler.polynomials import BaseStoredNASAPoly

if __name__ == '__main__':
    species = 'air'
    A = BaseStoredNASAPoly(species)
    print(A.cp)
