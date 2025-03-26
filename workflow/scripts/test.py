"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak
"""

import pytest
import msprime
from io_handler import MSerror
from msprime_simulation import (
    get_chromosome_bounds,
    create_recombination_map,
)

def test_get_chromosome_bounds():
    assert get_chromosome_bounds(1000) == (0, 1000)
    assert get_chromosome_bounds(500000) == (0, 500000)
    
    with pytest.raises(MSerror, match="Chromosome length cannot be negative."):
        get_chromosome_bounds(-100)

def test_create_recombination_map():
    chrom_length = 100000
    recombination_rate = 1e-8
    recomb_map = create_recombination_map(chrom_length, recombination_rate)
    
    assert isinstance(recomb_map, msprime.RateMap)
    assert recomb_map.num_intervals == 1