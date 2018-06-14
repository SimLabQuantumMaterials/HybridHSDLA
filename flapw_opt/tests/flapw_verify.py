#!/usr/bin/python
import numpy as np
import h5py as hp
import matplotlib.pyplot as plt

class HSDump:
    def __init__(self, fname, is_fleurdump=False):
        print "Opening %s with fleur=%s" % (fname, is_fleurdump)
        self.fname = fname
        self.is_fleurdump = is_fleurdump
        self.h5file = None
        self._reload_file()

    def _reload_file(self):
        if self.h5file:
            self.h5file.close()
        self.h5file = hp.File(self.fname)
        self.g_ham = self.h5file["matrix_hamilton"]
        self.g_olp = self.h5file["matrix_overlap"]
    def _get_h5_cx_set(self, dset):
        """Extract a complex dataset given in dset into a numpy array"""
        return np.squeeze(dset.value.view(np.complex))

    def read_hs(self, index):
        #s_raw = self.g_olp[str(index)].value
        #s = np.squeeze(s_raw.view(np.complex))
        s = self._get_h5_cx_set(self.g_olp[str(index)])
        h = self._get_h5_cx_set(self.g_ham[str(index)])
        if self.is_fleurdump:
            return pack_to_full(h), pack_to_full(s)
        else:
            return h.transpose(), s.transpose()
    def get_AB(self):
        if self.is_fleurdump:
            raise RuntimeError("Dataset does not contain A,B")

        A = self._get_h5_cx_set(self.h5file["A"])
        B = self._get_h5_cx_set(self.h5file["B"])
        return A, B

class Comparer:
    """d1 is the C++ dataset, d2 the one from Fleur"""
    def __init__(self, d1, d2):
        self.d1 = HSDump(d1)
        self.d2 = HSDump(d2, True)

    def reload_data(self):
        self.d1._reload_file()
        self.d2._reload_file()

    def __getitem__(self, i):
        return self.cmp(i)

    def cmp(self, i):
        h1, s1 = self.d1.read_hs(i)
        h2, s2 = self.d2.read_hs(i)
        return np.tril(h1 - h2), np.tril(s1 - s2)

    def print_error(self, imax=60):
        norm = np.linalg.norm
        self.reload_data()
        print "Index\t||hd||\t||sd||"
        for i in range(1,imax+1):
            hd, sd = self.cmp(i)
            print "%d:\t%g\t%g" % (i, norm(hd), norm(sd))


def pack_to_full(V):
    assert V.ndim == 1, "Need a vector, V.ndim = %d" % V.ndim
    l = V.size
    n = int(.5 * (np.sqrt(1+8*l) - 1))
    assert n*(n+1)/2 == l, "length does not yield a square matrix"
    M = np.zeros((n,n), V.dtype, order='F')
    row = col = 0
    for el in V:
        M[row, col] = el
        if (row == col):
            col = 0
            row += 1
        else:
            col += 1
    return M

def mplot(m):
    plt.imshow(np.abs(m), interpolation="nearest")
    plt.colorbar()
    return plt.show()
    

def get_stdcmp(exp_mine, exp_fleur):
    return Comparer(mpath.format(exp=exp_mine),
                    fpath.format(exp=exp_fleur))

cluster_home = "/home/mh447308"
cluster_work = "/work/mh447308"

vm_home = "/home/hrywniak"

# a little bit more sophisticated: check host, see if we have a better work dir guess
rootdirs = { 'rwth-aachen.de' : 
                    [cluster_work + '/th_cpp',
                     cluster_work + '/th_fleur'],
             'ubuntu' : # my VM
                    [vm_home + '/flapw_cpp/work',
                    vm_home + '/thesis/fleur/work']
            }
from platform import node
my_host = node().lower()
try:
    prefix = next(p for host, p in rootdirs.iteritems() if my_host.endswith(host))
except:
    raise LookupError("No prefix found for %s" % my_host)
    

mpath = prefix[0] + "/{exp}/flapw_out.h5"
fpath = prefix[1] + "/{exp}/fleur_dump.h5"

if __name__ == "__main__":
    import sys
    num_k = 1
    argc = len(sys.argv)
    if argc < 3:
        raise RuntimeError("Not enough arguments, need 2 (my_exp, fleur_exp)")
    exp_m = sys.argv[1]
    exp_f = sys.argv[2]
    if argc > 3:
        num_k = int(sys.argv[3])
    cmp = get_stdcmp(exp_m, exp_f)
    cmp.print_error(num_k)
    
