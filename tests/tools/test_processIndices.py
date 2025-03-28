import pytest
import datetime
from EUVpy.tools.processIndices import getCLSF107

def test_config():
    dateStart = "2000-01-01"
    dateEnd = "2000-01-31"
    times, F107, F107a, F107b = getCLSF107(dateStart, dateEnd, truncate=True)

    print(f'times: {times}, F107a: {F107a}, F107b: {F107b}')

    dateStartDatetime = datetime.datetime.strptime(dateStart, "%Y-%m-%d")
    dateEndDatetime = datetime.datetime.strptime(dateEnd, "%Y-%m-%d")
    numDays = (dateEndDatetime - dateStartDatetime).days
    assert(len(times) == numDays)
    assert(len(F107) == numDays)
    assert(len(F107a) == numDays)
    assert(len(F107b) == numDays)
    assert(F107a[0] == 161.14074074)