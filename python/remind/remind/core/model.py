#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the Model class."""

from __future__ import absolute_import

from cobra.core.object import Object
import cobra

import logging
import types

from collections import defaultdict
import pandas as pd
import optlang
import six
from optlang.symbolics import Basic, Zero
from six import iteritems
from optlang.exceptions import SolverError

from cobra.core.configuration import Configuration
from cobra.core.dictlist import DictList
from cobra.exceptions import SolverNotFound
from cobra.util.context import HistoryManager, resettable
from cobra.util.solver import (
    add_cons_vars_to_problem, assert_optimal, interface_to_str,
    remove_cons_vars_from_problem, set_objective, solvers)
from cobra.util.util import format_long_string

from pytfa.utils.str import camel2underscores
from pytfa.optim.utils import get_primal
from pytfa.core.model import timeit
from pytfa.optim.variables import GenericVariable

# from ..optim.utils import get_all_subclasses
from .solution import Solution


logger = logging.getLogger(__name__)
configuration = Configuration()


class Model(Object):
    """Class representation for a cobra model

    Parameters
    ----------
    id_or_model : Model, string
        Either an existing Model object in which case a new model object is
        instantiated with the same properties as the original model,
        or an identifier to associate with the model as a string.
    name : string
        Human readable name for the model

    Attributes
    ----------
    reactions : DictList
        A DictList where the key is the exchange identifier and the value a
        Reaction
    metabolites : DictList
        A DictList where the key is the metabolite identifier and the value a
        Metabolite
    species : DictList
        A DictList where the key is the species identifier and the value a
        Model
    groups : DictList
        A DictList where the key is the group identifier and the value a
        Group
    solution : Solution
        The last obtained solution from optimizing the model.

    """

    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model.
        """
        self.__dict__.update(state)
        for y in ['reactions', 'metabolites', 'species']:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __getstate__(self):
        """Get state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably.
        """
        odict = self.__dict__.copy()
        odict['_contexts'] = []
        return odict

    def __init__(self, id_or_model=None, name=None, sloppy=False):
        if isinstance(id_or_model, cobra.Model):
            Object.__init__(self, name=name)
            self.__setstate__(id_or_model.__dict__)
            if not hasattr(self, "name"):
                self.name = None
            self._solver = id_or_model.solver
        else:
            Object.__init__(self, id_or_model, name=name)

            self.species = DictList()
            self.reactions = DictList()  # A list of cobra.Reactions
            self.metabolites = DictList()  # A list of cobra.Metabolites
            self.groups = DictList()  # A list of cobra.Groups

            self._contexts = []

            # from cameo ...

            # Setting the optlang model
            interface = configuration.solver
            self._solver = interface.Model()
            self._solver.objective = interface.Objective(Zero)

            self._tolerance = None
            self.tolerance = configuration.tolerance
            
            self._cons_queue = list()
            self._var_queue = list()
    
            self._var_dict = dict()
            self._cons_dict = dict()
    
            self.sloppy=sloppy

    @property
    def solver(self):
        """Get or set the attached solver instance.

        The associated the solver object, which manages the interaction with
        the associated solver, e.g. glpk.

        This property is useful for accessing the optimization problem
        directly and to define additional non-metabolic constraints.

        Examples
        --------
        >>> import cobra.test
        >>> model = cobra.test.create_test_model("textbook")
        >>> new = model.problem.Constraint(model.objective.expression,
        >>> lb=0.99)
        >>> model.solver.add(new)
        """
        return self._solver

    @solver.setter
    @resettable
    def solver(self, value):
        not_valid_interface = SolverNotFound(
            '%s is not a valid solver interface. Pick from %s.' % (
                value, list(solvers)))
        if isinstance(value, six.string_types):
            try:
                interface = solvers[interface_to_str(value)]
            except KeyError:
                raise not_valid_interface
        elif isinstance(value, types.ModuleType) and hasattr(value, 'Model'):
            interface = value
        elif isinstance(value, optlang.interface.Model):
            interface = value.interface
        else:
            raise not_valid_interface

        # Do nothing if the solver did not change
        if self.problem == interface:
            return
        self._solver = interface.Model.clone(self._solver)

    @property
    def tolerance(self):
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value):
        solver_tolerances = self._solver.configuration.tolerances

        try:
            solver_tolerances.feasibility = value
        except AttributeError:
            logger.info("The current solver doesn't allow setting"
                        "feasibility tolerance.")

        try:
            solver_tolerances.optimality = value
        except AttributeError:
            logger.info("The current solver doesn't allow setting"
                        "optimality tolerance.")

        try:
            solver_tolerances.integrality = value
        except AttributeError:
            logger.info("The current solver doesn't allow setting"
                        "integrality tolerance.")

        self._tolerance = value
        
    def add_species(self, species):
        """Add a species to the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        species : Model
            A Model for the species to be added
        """
        # TODO
        raise NotImplementedError()

    def remove_species(self, species):
        """Remove a species from the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reactions : str
            An id for the model to be removed

        """
        # TODO
        raise NotImplementedError()

    def get_associated_groups(self, element):
        """Returns a list of groups that an element (reaction, metabolite, gene)
        is associated with.

        Parameters
        ----------
        element: `cobra.Reaction`, `cobra.Metabolite`, or `cobra.Gene`

        Returns
        -------
        list of `cobra.Group`
            All groups that the provided object is a member of
        """
        # check whether the element is associated with the model
        return [g for g in self.groups if element in g.members]

    def add_cons_vars(self, what, **kwargs):
        """Add constraints and variables to the model's mathematical problem.

        Useful for variables and constraints that can not be expressed with
        reactions and simple lower and upper bounds.

        Additions are reversed upon exit if the model itself is used as
        context.

        Parameters
        ----------
        what : list or tuple of optlang variables or constraints.
           The variables or constraints to add to the model. Must be of
           class `optlang.interface.Variable` or
           `optlang.interface.Constraint`.
        **kwargs : keyword arguments
           Passed to solver.add()
        """
        add_cons_vars_to_problem(self, what, **kwargs)

    def remove_cons_vars(self, what):
        """Remove variables and constraints from the model's mathematical
        problem.

        Remove variables and constraints that were added directly to the
        model's underlying mathematical problem. Removals are reversed
        upon exit if the model itself is used as context.

        Parameters
        ----------
        what : list or tuple of optlang variables or constraints.
           The variables or constraints to add to the model. Must be of
           class `optlang.interface.Variable` or
           `optlang.interface.Constraint`.
        """
        remove_cons_vars_from_problem(self, what)

    @property
    def problem(self):
        """The interface to the model's underlying mathematical problem.

        Solutions to cobra models are obtained by formulating a mathematical
        problem and solving it. Cobrapy uses the optlang package to
        accomplish that and with this property you can get access to the
        problem interface directly.

        Returns
        -------
        optlang.interface
            The problem interface that defines methods for interacting with
            the problem and associated solver directly.
        """
        return self.solver.interface

    @property
    def variables(self):
        """The mathematical variables in the cobra model.

        In a cobra model, most variables are reactions. However,
        for specific use cases, it may also be useful to have other types of
        variables. This property defines all variables currently associated
        with the model's problem.

        Returns
        -------
        optlang.container.Container
            A container with all associated variables.
        """
        return self.solver.variables

    @property
    def constraints(self):
        """The constraints in the cobra model.

        In a cobra model, most constraints are metabolites and their
        stoichiometries. However, for specific use cases, it may also be
        useful to have other types of constraints. This property defines all
        constraints currently associated with the model's problem.

        Returns
        -------
        optlang.container.Container
            A container with all associated constraints.
        """
        return self.solver.constraints

    @timeit
    def slim_optimize(self, error_value=float('nan'), message=None):
        """Optimize model without creating a solution object.

        Creating a full solution object implies fetching shadow prices and
        flux values for all reactions and metabolites from the solver
        object. This necessarily takes some time and in cases where only one
        or two values are of interest, it is recommended to instead use this
        function which does not create a solution object returning only the
        value of the objective. Note however that the `optimize()` function
        uses efficient means to fetch values so if you need fluxes/shadow
        prices for more than say 4 reactions/metabolites, then the total
        speed increase of `slim_optimize` versus `optimize` is  expected to
        be small or even negative depending on how you fetch the values
        after optimization.

        Parameters
        ----------
        error_value : float, None
           The value to return if optimization failed due to e.g.
           infeasibility. If None, raise `OptimizationError` if the
           optimization fails.
        message : string
           Error message to use if the model optimization did not succeed.

        Returns
        -------
        float
            The objective value.
        """
        self.solver.optimize()
        if self.solver.status == optlang.interface.OPTIMAL:
            return self.solver.objective.value
        elif error_value is not None:
            return error_value
        else:
            assert_optimal(self, message)

    def optimize(self, objective_sense=None, raise_error=False):
        """
        Optimize the model using flux balance analysis.

        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used.
        raise_error : bool
            If true, raise an OptimizationError if solver status is not
             optimal.

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.

        """
        if objective_sense:
            self.objective.direction = objective_sense
        
        original_direction = self.objective.direction
        self.objective.direction = \
            {"maximize": "max", "minimize": "min"}.get(
                objective_sense, original_direction)
        try:
            # self._hidden_optimize_call(kwargs)
            self.slim_optimize()
            solution = self.get_solution()
            self.solution = solution
            return solution
        except SolverError as SE:
            status = self.solver.status
            self.logger.error(SE)
            self.logger.warning('Solver status: {}'.format(status))
            raise (SE)

    def _repair(self, rebuild_index=True):
        """Update all indexes and pointers in a model

        Parameters
        ----------
        rebuild_index : bool
            rebuild the indices kept in reactions, metabolites and species
        rebuild_relationships : bool
             reset all associations between reactions, metabolites, model and
             then re-add them.
        """
        if rebuild_index:  # DictList indexes
            self.reactions._generate_index()
            self.metabolites._generate_index()
            self.species._generate_index()
            self.groups._generate_index()
        # if rebuild_relationships:
        #     for met in self.metabolites:
        #         met.species.clear()
        #         # met._reactions.clear()
        #     for rxn in self.reactions:
        #         rxn.species.clear()

        # point _model to self
        for l in (self.reactions, self.species, self.metabolites, self.groups):
            for e in l:
                e._model = self
    
    def repair(self):
        """
        Updates references to variables and constraints
        :return:
        """
        # self.add_cons_vars([x.constraint for x in self._cons_dict.values()])
        # self.add_cons_vars([x.variable for x in self._var_dict.values()])
        self._push_queue()
        self._repair()
        self.regenerate_constraints()
        self.regenerate_variables()

    @property
    def objective(self):
        """Get or set the solver objective

        Before introduction of the optlang based problems,
        this function returned the objective reactions as a list. With
        optlang, the objective is not limited a simple linear summation of
        individual reaction fluxes, making that return value ambiguous.
        Henceforth, use `cobra.util.solver.linear_reaction_coefficients` to
        get a dictionary of reactions with their linear coefficients (empty
        if there are none)

        The set value can be dictionary (reactions as keys, linear
        coefficients as values), string (reaction identifier), int (reaction
        index), Reaction or problem.Objective or sympy expression
        directly interpreted as objectives.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self.solver.objective

    @objective.setter
    def objective(self, value):
        if isinstance(value, Basic):
            value = self.problem.Objective(value, sloppy=False)
        if not isinstance(value, (dict, optlang.interface.Objective)):
            try:
                reactions = self.reactions.get_by_any(value)
            except KeyError:
                raise ValueError('invalid objective')
            value = {rxn: 1 for rxn in reactions}
        set_objective(self, value, additive=False)

    @property
    def objective_direction(self):
        """
        Get or set the objective direction.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when exiting the context.

        """
        return self.solver.objective.direction

    @objective_direction.setter
    @resettable
    def objective_direction(self, value):
        value = value.lower()
        if value.startswith("max"):
            self.solver.objective.direction = "max"
        elif value.startswith("min"):
            self.solver.objective.direction = "min"
        else:
            raise ValueError("Unknown objective direction '{}'.".format(value))

    def __enter__(self):
        """Record all future changes to the model, undoing them when a call to
        __exit__ is received"""

        # Create a new context and add it to the stack
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]

        return self

    def __exit__(self, type, value, traceback):
        """Pop the top context manager and trigger the undo functions"""
        context = self._contexts.pop()
        context.reset()

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Name</strong></td>
                <td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Number of metabolites</strong></td>
                <td>{num_metabolites}</td>
            </tr><tr>
                <td><strong>Number of reactions</strong></td>
                <td>{num_reactions}</td>
            </tr><tr>
                <td><strong>Number of groups</strong></td>
                <td>{num_groups}</td>
            </tr><tr>
                <td><strong>Objective expression</strong></td>
                <td>{objective}</td>
            </tr><tr>
                <td><strong>Compartments</strong></td>
                <td>{compartments}</td>
            </tr>
          </table>""".format(
            name=self.id,
            address='0x0%x' % id(self),
            num_metabolites=len(self.metabolites),
            num_reactions=len(self.reactions),
            num_groups=len(self.groups),
            objective=format_long_string(str(self.objective.expression), 100),
            compartments=", ".join(
                v if v else k for k, v in iteritems(self.compartments)
            ))

    def print_info(self):
        """
        Print information and counts for the cobra_model
        :return:
        """

        n_metabolites = len(self.metabolites)
        n_reactions = len(self.reactions)
        n_species = len(self.species)
        n_constraints = len(self.constraints)
        n_variables = len(self.variables)

        info = pd.DataFrame(columns=['value'])
        info.loc['name'] = self.name
        info.loc['description'] = self.description
        info.loc['num constraints'] = n_constraints
        info.loc['num variables'] = n_variables
        info.loc['num metabolites'] = n_metabolites
        info.loc['num reactions'] = n_reactions
        info.loc['num species'] = n_species
        info.index.name = 'key'

        print(info)

    def add_variable(self, kind, hook, queue=False, **kwargs):
        """ Add a new variable to a COBRApy cobra_model.

        :param kind:
        :param string,cobra.Reaction hook: Either a string representing the name
            of the variable to add to the cobra_model, or a reaction object if the
            kind allows it

        :returns: The created variable
        :rtype: optlang.interface.Variable

        """

        # Initialisation links to the cobra_model
        var = kind(hook,
                   # lb=lower_bound if lower_bound != float('-inf') else None,
                   # ub=upper_bound if upper_bound != float('inf') else None,
                   queue=queue,
                   **kwargs)

        self._var_dict[var.name] = var
        self.logger.debug('Added variable: {}'.format(var.name))
        # self.add_cons_vars(var.variable)

        return var

    def add_constraint(self, kind, hook, expr, queue=False,**kwargs):
        """ Add a new constraint to a COBRApy cobra_model

        :param kind:
        :param string,cobra.Reaction hook: Either a string representing the name
            of the variable to add to the cobra_model, or a reaction object if the
            kind allows it
        :param sympy.thermo.expr.Expr expr: The expression of the constraint

        :returns: The created constraint
        :rtype: optlang.interface.Constraint

        """

        if isinstance(expr, GenericVariable):
            # make sure we actually pass the optlang variable
            expr = expr.variable

        # Initialisation links to the cobra_model
        cons = kind(hook, expr, # problem = self.problem,
                    # lb=lower_bound if lower_bound != float('-inf') else None,
                    # ub=upper_bound if upper_bound != float('inf') else None,
                    queue=queue,
                    **kwargs)
        self._cons_dict[cons.name] = cons
        self.logger.debug('Added constraint: {}'.format(cons.name))
        # self.add_cons_vars(cons.constraint)

        return cons

    def _remove_associated_consvar(self, all_cons_subclasses, all_var_subclasses,
                                   collection):
        """
        Removes both the constraints and variables associated to an element,
        as long as it was used as a hook in the cons/var declaration.
        For example, upon removing a reaction, also removes its associated
        deltaG variables and coupling constraints
        """

        if not hasattr(collection, '__iter__'):
            collection = [collection]

        strfy = lambda x:x if isinstance(x, str) else x.id

        for cons_type in all_cons_subclasses:
            for element in collection:
                try:
                    cons = self._cons_kinds[cons_type.__name__].get_by_id(strfy(element))
                    self.remove_constraint(cons)
                except KeyError:
                    pass
        for var_type in all_var_subclasses:
            for element in collection:
                try:
                    var = self._var_kinds[var_type.__name__].get_by_id(strfy(element))
                    self.remove_variable(var)
                except KeyError:
                    pass


    def remove_variable(self, var):
        """
        Removes a variable

        :param var:
        :return:
        """
        # Get the pytfa var object if an optlang variable is passed
        if isinstance(var,optlang.Variable):
            var = self._var_dict[var.name]

        self._var_dict.pop(var.name)
        self.remove_cons_vars(var.variable)
        self.logger.debug('Removed variable {}'.format(var.name))

    def remove_constraint(self, cons):
        """
        Removes a constraint

        :param cons:
        :return:
        """
        # Get the pytfa var object if an optlang variable is passed
        if isinstance(cons,optlang.Constraint):
            cons = self._cons_dict[cons.name]

        self._cons_dict.pop(cons.name)
        self.remove_cons_vars(cons.constraint)
        self.logger.debug('Removed constraint {}'.format(cons.name))

    def _push_queue(self):
        """
        updates the constraints and variables of the model with what's in the
        queue
        :return:
        """

        self.add_cons_vars(self._var_queue, sloppy=self.sloppy)
        self.add_cons_vars(self._cons_queue, sloppy = self.sloppy)

        if len(self._var_queue) > 0:
            self.regenerate_variables()
        if len(self._cons_queue) > 0:
            self.regenerate_constraints()

        self._var_queue = list()
        self._cons_queue = list()


    def regenerate_variables(self):
        """
        Generates references to the cobra_model's constraints in self._var_dict
        as tab-searchable attributes of the thermo cobra_model
        :return:
        """

        # Let us not forget to remove fields that might be empty by now
        if hasattr(self, '_var_kinds'):
            for k in self._var_kinds:
                attrname = camel2underscores(k)
                try:
                    delattr(self, attrname)
                except AttributeError:
                    pass # The attribute may not have been set up yet

        _var_kinds = defaultdict(DictList)
        for k, v in self._var_dict.items():
            _var_kinds[v.__class__.__name__].append(v)

        for k in _var_kinds:
            attrname = camel2underscores(k)
            setattr(self, attrname, _var_kinds[k])

        self._var_kinds = _var_kinds

    def regenerate_constraints(self):
        """
        Generates references to the cobra_model's constraints in self._cons_dict
        as tab-searchable attributes of the thermo cobra_model
        :return:
        """

        # Let us not forget to remove fields that migh be empty by now
        if hasattr(self, '_cons_kinds'):
            for k in self._cons_kinds:
                attrname = camel2underscores(k)
                try:
                    delattr(self, attrname)
                except AttributeError:
                    pass # The attribute may not have been set up yet

        _cons_kinds = defaultdict(DictList)

        for k, v in self._cons_dict.items():
            _cons_kinds[v.__class__.__name__].append(v)

        for k in _cons_kinds:
            attrname = camel2underscores(k)
            setattr(self, attrname, _cons_kinds[k])

        self._cons_kinds = _cons_kinds

    def get_primal(self, vartype, index_by_reactions=False):
        """
        Returns the primal value of the cobra_model for variables of a given type

        :param index_by_reactions:
        :param vartype: Class of variable. Ex: pytfa.optim.variables.ThermoDisplacement
        :return:
        """
        return get_primal(self, vartype, index_by_reactions)

    def get_solution(self):
        """
        Overrides the cobra.thermo.solution method, to also get the supplementary
        variables we added to the cobra_model

        *   :code:`solution.fluxes` in `cobrapy` is a transformed version of the solver
            output, as it actually calculates the _net_ flux of each reaction by
            substracting the reverse variable value to the forward variable value.
            This should be used anytime one needs the actual flux value

        *   :code:`solution.raw` is a clear copy of the solver output. From there one
            can access the value at solution for all the variables of the problem.
            However, looking for a reaction ID in there will only give the
            _forward_ flux. This should be used for any other variable than fluxes.

        *   :code:`solution.values` yields variables multiplied by their scaling factor
            (1 by default). Useful if you operated scaling on your equations for
            numerical reasons. This does _not_ include fluxes

        :return:
        """
        objective_value = self.solver.objective.value
        status = self.solver.status
        variables = pd.Series(data=self.solver.primal_values)

        solution = Solution(objective_value=objective_value, status=status,)

        self.solution = solution

        self.solution.raw = variables

        self.\
            solution.values = pd.DataFrame.from_dict({k:v.unscaled
                                                 for k,v in self._var_dict.items()},
                                                 orient = 'index')

        return solution

    def get_constraints_of_type(self, constraint_type):
        """
        Convenience function that takes as input a constraint class and returns
        all its instances within the cobra_model

        :param constraint_type:
        :return:
        """
        if isinstance(constraint_type,str):
            constraint_key = constraint_type
        else:
            #it is a class
            constraint_key = constraint_type.__name__
        return self._cons_kinds[constraint_key]

    def get_variables_of_type(self, variable_type):
        """
        Convenience function that takes as input a variable class and returns
        all its instances within the cobra_model

        :param variable_type:
        :return:
        """
        if isinstance(variable_type,str):
            variable_key = variable_type
        else:
            #it is a class
            variable_key = variable_type.__name__
        return self._var_kinds[variable_key]
