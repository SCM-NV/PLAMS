#!/usr/bin/env amspython
# coding: utf-8

import pytest

from scm.plams.tools.postprocess_results import moving_average


def test_moving_average():
    avg_x, avg_y = moving_average([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0], window=2)

    assert avg_x.tolist() == pytest.approx([1.5, 2.5, 3.5])
    assert avg_y.tolist() == pytest.approx([3.0, 5.0, 7.0])

    avg_x, avg_y = moving_average([1.0, 2.0], [3.0, 5.0], window=0)

    assert avg_x.tolist() == pytest.approx([1.0, 2.0])
    assert avg_y.tolist() == pytest.approx([3.0, 5.0])
