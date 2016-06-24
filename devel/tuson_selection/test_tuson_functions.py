import simuOpt
simuOpt.setOptions(quiet=True, numThreads=4)
import simuPOP as sim
import os, numpy as np, pandas as pd, collections as col
from saegus import analyze, simulate, parameters

import pytest

@pytest.fixture(scope="session")
def tuson_pop():
    tuson = sim.loadPopulation('tuson.pop')
    return tuson

def test_info_fields(tuson_pop):
    pop = tuson_pop
    assert 'primary' in pop.infoFields(), "Population lacks the primary infoField. Cannot assign structured mating."