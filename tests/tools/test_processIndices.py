import pytest
import datetime
from EUVpy.tools.processIndices import getCLSF107
#-----------------------------------------------------------------------------------------------------------------------

# Helper function:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    # https://stackoverflow.com/questions/5595425/how-to-compare-floats-for-almost-equality-in-python
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

#-----------------------------------------------------------------------------------------------------------------------
def test_config():
    dateStart = "2000-01-01"
    dateEnd = "2000-01-31"
    times, F107, F107a, F107b = getCLSF107(dateStart, dateEnd, truncate=False)

    print(f'times: {times}, F107a: {F107a}, F107b: {F107b}')

    dateStartDatetime = datetime.datetime.strptime(dateStart, "%Y-%m-%d")
    dateEndDatetime = datetime.datetime.strptime(dateEnd, "%Y-%m-%d")
    numDays = (dateEndDatetime - dateStartDatetime).days
    assert(len(times) == numDays)
    assert(len(F107) == numDays)
    assert(len(F107a) == numDays)
    assert(len(F107b) == numDays)
    # assert(F107a[0] == 161.14074074)
    assert(isclose(F107a[0], 161.14074074, rel_tol=1e-7))