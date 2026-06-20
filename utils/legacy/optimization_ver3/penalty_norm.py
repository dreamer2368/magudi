import numpy as np

__all__ = ['PenaltyBase','HuberNorm','penaltySwitcher']

class PenaltyBase:
    def norm(self, L2square):
        return 0.5 * L2square

    def gradientFactor(self, L2square):
        return 1.0

class HuberNorm(PenaltyBase):
    def __init__(self, eps=1.0e-6):
        self.eps = eps

    def norm(self, L2square):
        L2norm = np.sqrt( L2square )
        return ( (L2norm - 0.5*self.eps) if (L2norm>=self.eps) else (0.5 * L2square / self.eps) )

    def gradientFactor(self, L2square):
        L2norm = np.sqrt( L2square )
        return ( L2norm if (L2norm>=self.eps) else self.eps )

penaltySwitcher = {
    'base': PenaltyBase(),
    'huber': HuberNorm(),
}
