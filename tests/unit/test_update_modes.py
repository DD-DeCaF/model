import logging
from tempfile import mkdtemp
import os

from cobra.io import read_sbml_model

from model.adapter import get_unique_metabolite
from model.update_models import update_local_models

logging.disable(logging.CRITICAL)

def test_update_models():
    tempdir = mkdtemp()
    update_local_models('e_coli_core', tempdir)
    model = read_sbml_model(os.path.join(tempdir, 'e_coli_core.sbml.gz'))
    glucose = get_unique_metabolite(model, 'CHEBI:42758')
    assert glucose.id == 'glc__D_e'
    assert glucose.annotation['bigg.metabolite'] == 'glc__D'
