# DEQscaler
"""
Tool to rescale systems of ordinary differential equations of the type dy /dt = f(t,y) with initial conditions y(t0)=y0 and parameters.
1. A DEQ is formulated in a sympy representation.
2. The sympy representation is used to apply the parameters and thence create f(t,y) only.
3. The equation dy/dt=f(t,y) is solved numerically according to given precision and method by scipy.integrate.solve_ivp
4. The absolut maximum values of the numerical solution are used to rescale the problem to the domain [-1,1] (analog computer implementation).
"""

import numpy as np
import sympy
from scipy import integrate

"""
A class to represent and rescale a system of ordinary differential equations of the type dy/dt = f(t,y) with initial conditions y(t0)=y0 and parameters.
The rescaling is done by determining the maxima of the trajectory components of y using the solve_ivp function from scipy.

args_solve_ivp_ini_val: List of essential arguments t_span and y0 for scipy.integrate.solve_ivp. Cf. for further detail: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
t:                      sympy.Symbol('t') the time variable.
dydt:                   List of sympy symbols representing the components of y. 
f_t_y:                  List of sympy expressions representing f(t,y) including the sympy symbols representing the parameters. Note that t may not be part of f(t,y) at this time.
diff_eq_parameters:     Dictionary with the parameter sympy symbols as keys and their values.
max_scale_factor:       Optional factor to adapt the maxima used for rescaling.
kwargs_solve_ivp:       Optional arguments that can be passed to scipy.integrate.solve_ivp. Cf. link above for details.    
"""
class DEQscaler:

    def __init__(self, args_solve_ivp_ini_val, t, dydt, f_t_y, diff_eq_parameters, max_scale_factor=1.0,
                 kwargs_solve_ivp={}):
        self.args_solve_ivp_ini_val = args_solve_ivp_ini_val
        self.t = t
        self.dydt = dydt
        self.f_t_y = f_t_y
        self.diff_eq_parameters = diff_eq_parameters
        self.max_scale_factor = max_scale_factor
        self.kwargs_solve_ivp = kwargs_solve_ivp

        self.num_f_t_y = []
        self.num_sol = []
        self.maxima = {}

    """
    Display the properties of the created object as text in the console.
    """
    def show_eqn(self):
        print("\nInitial values: Starting value y(t0)=y0; integrate numerically wrt. t until tf is reached.")
        print("y0       = " + str(self.args_solve_ivp_ini_val[1]))
        print("(t0, tf) = " + str(self.args_solve_ivp_ini_val[0]))
        if self.max_scale_factor != 1.0:
            print("max_scale_factor = " + str(self.max_scale_factor))

        if self.kwargs_solve_ivp != {}:
            print("\nAdditional kwargs passed to scipy.integrate.solve_ivp:")
            print(self.kwargs_solve_ivp)

        print("")

        print("Parameters:")
        for label, val in self.diff_eq_parameters.items():
            print("{:} = {:}".format(str(label), val))

        print("")

        print("System of Equations: [Note the order is alphabetical in the sympy-symbols present.]")
        for pos, ele in enumerate(self.dydt):
            print("{:} = {:}".format(str(ele), str(self.f_t_y[pos])))

        print("")

        print("System of Equations (for numerical calculation):")
        """Create tuple to substitute parameters with their values."""
        substitution_parameters = [(ele[0], ele[1]) for ele in self.diff_eq_parameters.items()]
        """Use the substitution to display the numerical values."""
        for pos, ele in enumerate(self.dydt):
            print("{:} = {:}".format(str(ele), str(self.f_t_y[pos].subs(substitution_parameters))))

    """
    Use sympy.lambdify to create numeric scipy functions with substituted values for the parameters.
    """
    def create_num_f_t_y(self):

        """
        Create tuple of the variables used to solve later on.
        """
        calc_variables = tuple([self.t] + self.dydt)

        """
        Iterate over the f(t,y) functions, i.e. proceed elementwise.
        For every function use the .subs() to replace the parameters with their numeric values.        
        Use the tuple calc_variables to specify on which variables the numeric functions depend.        
        Ensure a scipy object is created using the modules kwarg.
        """
        self.num_f_t_y = [
            sympy.lambdify([calc_variables], ele.subs(self.diff_eq_parameters), modules=["scipy", "numpy"]) for ele in
            self.f_t_y]

    """
    The scipy.integrate.solve_ivp requires a callable f(t,y) function fun.
    This is implemented here using the num_f_t_y.
    """
    def derivative(self, t, state):
        """
        Create List of step values including t. [The t might not be used according to the specific equation, but is required for solve_ivp.]
        """
        var_step_vals = np.concatenate(([t], state))

        """
        Evaluate num_f_ty using the step values.
        """
        calc_eqns = [ele(var_step_vals) for ele in self.num_f_t_y]

        return calc_eqns

    """
    Perform the numeric integration to solve the system using the args and potential kwargs specified.
    """
    def solve_numerically(self):
        """
        Create the lambdify numeric functions to be used
        """
        self.create_num_f_t_y()

        """
        Use the scipy.integrate.solve_ivp function.
        """
        self.num_sol = integrate.solve_ivp(self.derivative, *self.args_solve_ivp_ini_val, **self.kwargs_solve_ivp)

    """
    Determine the absolut maxima of the trajectories using the numerical solution.
    """
    def determine_max(self):
        """
        Calculate the numerical solution.
        """
        self.solve_numerically()

        """
        Determine the absolut maxima.
        Potential adaption: Use the numpy ceiling function np.ceil to round to 3 significant figures.
        by replacint "np.absolute(ele).max()" with "np.ceil(np.absolute(ele).max()*1000)/1000"
        """
        self.maxima = {self.dydt[pos]: np.absolute(ele).max() for pos, ele in enumerate(self.num_sol.y)}

    """
    Create arguments to instanciate the new, rescaled DEQ system 
    """
    def create_rescaled_Diff_Eq(self, maxima_items=None):
        """
        If no individual maxima are provided calculate maxima by solving numerically and use these.
        """
        if maxima_items is None:
            self.determine_max()
            maxima_items = self.maxima

        """
        Make maxima_items iterable.
        """
        maxima_items = maxima_items.items()

        """
        Adapt the maxima by multipling with the max_scale_factor.
        """
        round_up_maxima = {key: value * self.max_scale_factor for key, value in maxima_items}

        """
        Rescale the initial values, i.e. y0 in the args_solve_ivp_ini_val
        """
        rescaled_args_solve_ivp_ini_val = [self.args_solve_ivp_ini_val[0],
                                           [self.args_solve_ivp_ini_val[1][pos] * 1 / round_up_maxima[ele] for pos, ele
                                            in enumerate(self.dydt)]]

        """
        Create the tuples used to substitute the symbolic variables with the product of the rescaling factor and the symbolic variable.
        """
        substitution_list_tuples = [(ele[0], ele[0] * ele[1]) for ele in round_up_maxima.items()]

        """
        Create the factor list used to scale each expression in f_t_y as a whole.
        """
        rescale_eqns_factors = [1 / ele[1] for ele in round_up_maxima.items()]

        """
        Create the rescaled equation rescaled_f_t_y.
        """
        rescaled_f_t_y = [rescale_eqns_factors[pos] * ele.subs(substitution_list_tuples) for pos, ele in
                          enumerate(self.f_t_y)]

        return [rescaled_args_solve_ivp_ini_val, self.t, self.dydt, rescaled_f_t_y, self.diff_eq_parameters,
                self.max_scale_factor, self.kwargs_solve_ivp]
