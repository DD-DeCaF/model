import logging
import pytest

from model.adapter import NoIDMapping, full_genotype, get_unique_metabolite
from model.storage import find_in_memory

logging.disable(logging.CRITICAL)


@pytest.mark.asyncio
async def test_existing_metabolite():
    ecoli = find_in_memory('iJO1366')
    assert get_unique_metabolite(ecoli, 'chebi:17790') == get_unique_metabolite(
        ecoli, 'meoh', db_name='bigg.metabolite')
    assert get_unique_metabolite(ecoli, 'succ', db_name='bigg.metabolite').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        await get_unique_metabolite(ecoli, 'wrong_id')
