#Test polyfit function used for linear calibration by showiung that the midpoint height (avg of y-values) is equal to calculated midpoint (avg of x-values * slope + intercept)

from nose.tools import assert_almost_equal
import numpy as np
import pytest

def test_polyfit():
    x = [45, 513]
    y = [15, 38]
    slope, intercept = np.polyfit(x, y, 1)
    midpoint_x = (45+513)/2
    midpoint_y = (15+38)/2
    assert_almost_equal(midpoint_x*slope + intercept, midpoint_y)
