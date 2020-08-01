import scipy as sp
import ot

def gw(xs, xt, n_samples):
  # compute distance kernels 
  C1 = sp.spatial.distance.cdist(xs, xs)
  C2 = sp.spatial.distance.cdist(xt, xt)
  # normalize distance kernels
  C1 /= C1.max()
  C2 /= C2.max()
  # Compute Gromov-Wasserstein distance
  p = ot.unif(n_samples)
  q = ot.unif(n_samples)
  return ot.gromov.gromov_wasserstein(C1, C2, p, q, 'square_loss', verbose=False, log=True)
