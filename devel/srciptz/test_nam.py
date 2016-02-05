import simuOpt
simuOpt.setOptions(quiet=True, numThreads=4)
import simuPOP as sim
import pytest

@pytest.fixture
def load_nam():
    nam = sim.loadPopulation('nam_prefounders.pop')
    return nam

def test_loaded_nam(load_nam):
    size = load_nam.popSize()
    assert size == 1 "Population size should be 26"

def test_
    
    

    