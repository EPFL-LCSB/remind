# -*- coding: utf-8 -*-

"""Provide unified interfaces to optimization solutions."""

from __future__ import absolute_import

import logging
from builtins import object, super

from optlang.interface import OPTIMAL
from pandas import option_context

LOGGER = logging.getLogger(__name__)


class Solution(object):
    """
    A unified interface to a `cobra.Model` optimization solution.

    Notes
    -----
    Solution is meant to be constructed by `get_solution` please look at that
    function to fully understand the `Solution` class.

    Attributes
    ----------
    objective_value : float
        The (optimal) value for the objective function.
    status : str
        The solver status related to the solution.
    reduced_costs : pandas.Series
        Contains reaction reduced costs (dual values of variables).
    shadow_prices : pandas.Series
        Contains metabolite shadow prices (dual values of constraints).
    """

    def __init__(self, objective_value, status, reduced_costs=None,
                 shadow_prices=None, **kwargs):
        """
        Initialize a `Solution` from its components.

        Parameters
        ----------
        objective_value : float
            The (optimal) value for the objective function.
        status : str
            The solver status related to the solution.
        reduced_costs : pandas.Series
            Contains reaction reduced costs (dual values of variables).
        shadow_prices : pandas.Series
            Contains metabolite shadow prices (dual values of constraints).
        """
        super(Solution, self).__init__(**kwargs)
        self.objective_value = objective_value
        self.status = status
        self.reduced_costs = reduced_costs
        self.shadow_prices = shadow_prices

    def __repr__(self):
        """String representation of the solution instance."""
        if self.status != OPTIMAL:
            return "<Solution {0:s} at 0x{1:x}>".format(self.status, id(self))
        return "<Solution {0:.3f} at 0x{1:x}>".format(self.objective_value,
                                                      id(self))

    def _repr_html_(self):
        if self.status == OPTIMAL:
            with option_context('display.max_rows', 10):
                html = ('<strong><em>Optimal</em> solution with objective '
                        'value {:.3f}</strong><br>{}'
                        .format(self.objective_value,
                                self.to_frame()._repr_html_()))
        else:
            html = '<strong><em>{}</em> solution</strong>'.format(self.status)
        return html
