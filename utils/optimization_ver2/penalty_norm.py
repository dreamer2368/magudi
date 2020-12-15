import numpy as np

__all__ = ['penaltyBase','HuberNorm','penaltySwitcher']

class penaltyBase:
    def norm(L2square):
        return 0.5 * L2square

    def gradientFactor(L2square):
        return 1.0

class HuberNorm(penaltyBase):
    def __init__(self, eps=1.0e-6):
        self.eps = eps

    def norm(L2square):
        L2norm = np.sqrt( L2square )
        return ( (L2norm - 0.5*self.eps) if (L2norm>=self.eps) else (0.5 * L2square / self.eps) )

    def gradientFactor(L2square):
        L2norm = np.sqrt( L2square )
        return ( L2norm if (L2norm>=self.eps) else self.eps )

penaltySwitcher = {
    'base': penaltyBase(),
    'Huber': HuberNorm(),
}
