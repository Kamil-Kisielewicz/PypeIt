# Module to run tests on ararclines
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pytest

from pypeit import ararclines


def test_load_linelist():
    """ Touches nearly all the methods in ararclines
    Returns
    -------
    """
    # Init
    # Load
    alist = ararclines.load_arcline_list(['CuI','ArI','NeI'], '600/7500', 'keck_lris_red',
                                         modify_parse_dict=dict(NeI={'min_wave': 3000.},
                                                                ArI={'min_intensity': 399.}))
    # Min NeI
    NeI = alist['Ion'] == 'NeI'
    np.testing.assert_allclose(np.min(alist['wave'][NeI]), 3455.1837999999998)

