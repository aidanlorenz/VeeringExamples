#!/usr/bin/env python
# coding: utf-8

# In[2]:


# IF CHANGES ARE MADE IN THIS FILE, IT WILL NEED TO BE REDOWNLOADED AS A .PY (or maybe .ipynb) FILE
# THEN CHANGE THE EXTENSION TO .SAGE AND MOVED TO THE GENERAL FOLDER FOR THE OTHER FILES TO LOAD IT
# PROPERLY.

from veering import taut_polytope
from veering import taut
from veering import taut_polynomial
from veering import transverse_taut
from veering import carried_surface
import snappy

from sage.plot.contour_plot import ContourPlot

import numpy as np

import time

m = mathematica


# In[3]:


sig = 'kLLLPLQkcefegijjiijiieldllxtxa_0112222011'
# SOME INTERESTING EXAMPLES ON WHICH TO TEST

# Homology dimension 1 example:
#sig = 'cPcbbbiht_12'
# Homology dimension 2 example: (Hironaka example)
#sig = 'eLMkbcddddedde_2100'
# Homology dimension 3 example:
#sig = 'ivLLQQccdgffhghhvvaaaaavv_01111220'
# Homology dimension 4 example:
#sig = 'ivvPQQcfhghgfghfaaaaaaaaa_01122000'
# Homology with torsion example (dimension 3):
#sig = 'jvLPzQQedfgihghiivkqaaacvvc_102021211'

#sig = 'gvLQQcdeffeffffaafa_201102' # KKT Example

# Working Examples
#sig = 'fLLQcbecdeepuwsua_20102'
#sig = 'fLLQcbeddeehhbghh_01110'
#sig = 'jLLAvQQbcdehiihhitsjaeqmqut_201021210'

#sig = 'hLLMMkaedfdgggjxaxjxqh_2002110' # Example of a dimension 2 homology with only 1 boundary component.


# In[4]:


def get_hom_vars(sig):
    '''Auxiliary function which creates a tuple
    of variables of length given by the dimension
    of the homology for the 3-manifold given by
    the sig.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
        
    Return type:
        tuple(sage.symbolic.expression.Expression)
        
    Dependencies: None 
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    return tuple(var('a' + str(i)) for i in range(tri.homology().rank()))


# In[5]:


get_hom_vars(sig)


# In[6]:


def get_extreme_surfs(sig): # NOT SURE IF I ACTUALLY NEED THIS FUNCTION FOR ANYTHING...
    '''DEPRECATED Gets the branch equation 
    (non-negative 2-chain) form of the surfaces 
    which span the fibered cone.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list(sage.modules.vector_integer_dense.Vector_integer_dense)
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    extreme_rays = taut_polytope.cone_in_homology(sig) # rays (in homology) spanning the fibered cone.
    rays = taut_polytope.taut_rays(sig) # These are rays spanning the cone of carried surfaces (not sure exactly what these are tbh...)
    P = taut_polytope.projection_to_homology(sig) # This is the matrix to multiply by to get branched surface tuple to a second homology class.
    extreme_surfs = []
    for ray in rays: # this for loop stores the surfaces (in terms of their branch equations) that define the boundary of the cone
        if P*ray in extreme_rays:
            temp = 0
            for surf in extreme_surfs:
                if P*surf != P*ray:
                    temp = temp + 1
            if temp == len(extreme_surfs):
                extreme_surfs.append(ray)
            
    if tri.homology().rank() > 1 and [tuple(extreme_rays[0]), tuple(extreme_rays[1])] != [tuple(P*extreme_surfs[0]), tuple(P*extreme_surfs[1])]:
        extreme_surfs.reverse()
            
    return extreme_surfs


# In[7]:


get_extreme_surfs(sig)


# In[8]:


def get_fibered_cone(sig):
    '''Builds the fibered cone in homology
    of the given veering triangulation.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.geometry.cone.ConvexRationalPolyhedralCone
    '''
    extreme_rays = taut_polytope.cone_in_homology(sig)
    return Cone(extreme_rays)


# In[9]:


get_fibered_cone(sig)


# In[10]:


def get_fibered_cone_plot(sig):
    '''Plots the fibered cone from
    get_fibered_cone. It will only work
    if the homology has dimension 2 or 3.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.plot.graphics.Graphics
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    rank = tri.homology(1).rank()
    if rank == 2 or rank == 3:
        return get_fibered_cone(sig).plot(mode='box', xmin = -10, xmax = 10, ymin = -10, ymax = 10, zmin = -10, zmax = 10, ray_thickness=2, point_size=10)
    
    return 'This function only works for manifolds with H_2(M,dM) (equivalently H_1(M)) dimension 2 or 3.'


# In[11]:


get_fibered_cone_plot(sig)


# In[12]:


def get_spec(sig, method = 'fox'):
    '''Auxiliary function for the get_spec_at_e function.
    Returns the specialization of the taut polynomial at an
    arbitrary point in homology. Also used in 
    the get_dila_mathematica function.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        sage.symbolic.expression.Expression
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    tri, angle = taut.isosig_to_tri_angle(sig)
    #rank = tri.homology(1).rank()
    hom_vars = get_hom_vars(sig)
    if method == 'original':
        poly = taut_polynomial.taut_polynomial(tri, angle)
    elif method == 'fox':
        poly = taut_polynomial.taut_polynomial_via_fox_calculus(tri, angle)
    elif method == 'tree':
        poly = taut_polynomial.taut_polynomial_via_tree(tri, angle)
    unspec_monos = poly.monomials()
    monos = [] # the monomials in the specialization
    for i in range(poly.number_of_terms()):
        monos.append(unspec_monos[i].lc()*x^(sum([poly.exponents()[i][j]*hom_vars[j] for j in range(len(hom_vars))])))
        
    spec = 0
    for i in range(len(monos)):
        spec = spec + monos[i]*poly.monomial_coefficient(poly.monomials()[i])
        
    return spec


# In[13]:


get_spec(sig)


# In[14]:


def get_spec_at_e(sig, method = 'fox'):
    '''Returns the specialization (at an arb. point)
    of the taut polynomial with e plugged in. This
    amounts to substituting x = e in the output of
    get_spec. Used in the get_fibered_cone_and_levset_plot
    function as well as the get_minimal_direction function.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        sage.symbolic.expression.Expression
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    return get_spec(sig, method)(x = e)


# In[15]:


get_spec_at_e(sig)


# In[16]:


def get_fibered_cone_and_levset_plot(sig, method = 'fox'):
    '''Plots both the fibered cone and the level set
    where the dilatation is equal to e (so log of the
    dilatation is 1). This is significant because on
    this level set, the normalized dilatation is just
    e^{Thurston norm} (or just equal to the Thurston
    norm if using the log version of normalized dilatation.)
    Only works if the homology is dimension 2 or 3.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        sage.plot.graphics.Graphics
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    tri, angle = taut.isosig_to_tri_angle(sig)
    rank = tri.homology(1).rank()
    if rank == 3:
        levset = implicit_plot3d(get_spec_at_e(sig, method) == 0, (-10,10), (-10,10), (-10,10), cmap=['green'])
        Plot = levset + get_fibered_cone_plot(sig)
        return Plot
    
    if rank == 2:
        levset = implicit_plot(get_spec_at_e(sig, method) == 0, (-10,10), (-10,10), cmap=['green'])
        Plot = levset + get_fibered_cone_plot(sig)
        return Plot
    
    return 'This function only works for manifolds with H_2(M,dM) (equivalently H_1(M)) dimension 2 or 3.'


# In[17]:


get_fibered_cone_and_levset_plot(sig)


# In[18]:


def get_C2M_to_C1M(sig):
    '''Generates the standard boundary map
    on 2-chains, with orientations agreeing
    with the veering coorientations.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.matrix.matrix_integer_dense.Matrix_integer_dense
    '''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    tet_vert_coorientations = transverse_taut.is_transverse_taut(tri, angle, return_type = "face_coorientations") # This gives a list telling us whether the RHR orientations of the regina numbering agrees or disagrees with the veering coorientation. (See Veering repo for details)
    
    reg_mat = tri.boundaryMap(2) # This is the Regina generated boundary map, but the signs are according to traversing the triangles with the RHR wrt the Regina vertex ordering.

    nrows = reg_mat.rows()
    ncols = reg_mat.columns()
    
    # We first make a matrix
    mat = [[0 for i in range(ncols)] for j in range(nrows)]
    for i in range(nrows):
        for j in range(ncols):
            mat[i][j] = int(str(reg_mat.entry(i,j)))
        
    # And we change the signs to agree with the veering coorientation (multiply entire columns by the discrepancy factors +-1 recorded in tet_vert_coorientations)
    for i in range(ncols):
        for j in range(nrows):
            mat[j][i] = tet_vert_coorientations[i]*mat[j][i]
    
    return matrix(mat)


# In[19]:


get_C2M_to_C1M(sig)


# In[20]:


def get_2_chain_from_point_in_homology(sig):
    '''Produces a 2-chain which representing an
    arbitrary homology class, making sure that the
    obtained 2-chain is not only a preimage under
    the projection from 2-chains to homology, but
    also that it represents a homology class; i.e.
    it is in the kernel of the boundary map.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.modules.free_module.FreeModule_ambient_field_with_category.element_class
    '''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    hom_vars = get_hom_vars(sig)
    #print(hom_vars) # For testing
    
    P = taut_polytope.projection_to_homology(sig)
    #print(P) # For testing
    num_triangles = tri.countTriangles() # This is the number of variables we need to store for our 2-chain.
    
    C2M_to_C1M = get_C2M_to_C1M(sig)
    #chain_vars = var(','.join('a%s'%i for i in range(num_triangles))) # The 2-chain variables
    chain_vars = tuple(var('x' + str(i)) for i in range(num_triangles))
    preimage_in_kernel = vector(chain_vars)
    #print(preimage_in_kernel) # For testing
    
    eqns = []
    for term in list(C2M_to_C1M*preimage_in_kernel):
        eqns.append(term == 0) # we want our 2-chain to be in the kernel of the boundary map so that it actually represents a homology class
    #print(eqns) # For testing
    for i in range(len(hom_vars)):
        eqns.append(list(P*preimage_in_kernel)[i] == hom_vars[i]) # and we want its image in homology to be our homology class
    #print(eqns) # For testing
    solution = solve(eqns, chain_vars)
    
    preimage = vector([solution[0][i].rhs() for i in range(len(solution[0]))])
    
    return preimage


# In[21]:


get_2_chain_from_point_in_homology(sig)


# In[22]:


def get_Thurston_norm(sig):
    '''Gets the formula for the Thurston norm
    of our 3-manifold in the fibered face
    corresponding to the given veering triangulation
    for an arbitrary homology class.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.symbolic.expression.Expression
    '''
            
    surf = get_2_chain_from_point_in_homology(sig)
    # somehow the free variables from the get_2_chain_from_point_in_homology always all cancel so aren't left in the final formula
    prenorm = 0
    for tri in surf:
        prenorm = prenorm + tri
        
    return prenorm/2


# In[23]:


get_Thurston_norm(sig)


# In[24]:


def get_spec_with_subs(sig, point, method = 'fox'):
    '''Returns the specialization of the taut
    polynomial for a given point in homology.
    This amounts to replacing the variables in the
    output of get_spec(sig, method) with the numeric values
    in point. Used in the get_dila_mathematica
    function.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    point: tuple
        The point in H_2(M, dM; Z) at which you want to
        get the specialization.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)

    
    Return type:
        sage.symbolic.expression.Expression
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
    
    keys = get_hom_vars(sig)
    
    if len(keys) != len(point):
        return 'This point is not of the correct dimension.'
        
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitute = {} # Need to make a dictionary to keep track of which variables we are substituting for which variable, since the get_spec function returns the specialization with variables and we need to plug in our point.
    for key, coord in zip(keys, point):
        substitute[key] = coord
    
    poly = get_spec(sig, method).subs(substitute)
    
    return poly


# In[25]:


get_spec_with_subs(sig, (-1,2.9,1))


# In[26]:


def get_dila_mathematica(sig, point, method = 'fox'): # LOOK CAREFULLY INTO THIS FUNCTION...ALL CAPS COMMENT
    '''Returns the dilatation of the monodromy
    associated to point in the manifold
    specified by sig. Used in get_minimal_direction
    and get_min_norm_dila_log_mathematica and
    get_norm_dila_log_mathematica.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    point: tuple
        The point in H_2(M, dM; Z) at which you want to
        get the specialization.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        sage.symbolic.expression.Expression (really just a number though)
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
    
    keys = get_hom_vars(sig)
    
    if len(keys) != len(point):
        return 'This point is not of the correct dimension.'
        
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    #substitute = {} # Need to make a dictionary to keep track of which variables we are substituting for which variable, since the get_spec function returns the specialization with variables and we need to plug in our point.
    #for key, coord in zip(keys, point):
    #    substitute[key] = coord
    
    #poly = str(get_spec(sig).subs(substitute)) # substitute our point into the specialization with the variables
    poly = str(get_spec_with_subs(sig, point, method))
    #print(poly) # For testing
    poly = poly.replace("sqrt(x)", "Sqrt[x]") # NOT SURE WHEN THIS IS A PROBLEM...MAYBE THIS IS SOMETHING THAT NEEDS TO BE LOOKED INTO MORE GENERALLY?
    #print(poly) # For testing
    roots = m('roots = NSolve[' + poly + ' == 0, x]') # use mathematica to find the roots
    #print(roots) # For testing
    root_nums = []
    for i in range(len(roots)):
        root_nums.append(m('x/.roots[[' + str(i+1) + ']]')) # get the roots in a list, just of the values
        
    root_nums = [abs(root.sage()) for root in root_nums] # change them to sage objects and take their absolute values
    return max(root_nums) # and return the largest


# In[27]:


get_dila_mathematica(sig, (-1,2.9,1))


# In[28]:


def get_hom_dila_mathematica(sig, point): # WORK ON THIS
    '''DEPRECATED. This function takes in a point in
    homology and returns the homological
    dilatation of that point. Useful for seeing
    whether a given monodromy is orientable or
    not.'''
    
    keys = get_hom_vars(sig)
    
    if len(keys) != len(point):
        return 'This point is not of the correct dimension.'
        
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitute = {} # Need to make a dictionary to keep track of which variables we are substituting for which variable, since the get_spec function returns the specialization with variables and we need to plug in our point.
    for key, coord in zip(keys, point):
        substitute[key] = coord
    
    poly = str(snappy.Manifold(tri).alexander_polynomial().subs(substitute)) # NEEDS WORK HERE
    print(poly) # For testing
    poly = poly.replace("sqrt(x)", "Sqrt[x]") # NOT SURE WHEN THIS IS A PROBLEM...MAYBE THIS IS SOMETHING THAT NEEDS TO BE LOOKED INTO MORE GENERALLY?
    #print(poly) # For testing
    roots = m('roots = NSolve[' + poly + ' == 0, x]') # use mathematica to find the roots
    #print(roots) # For testing
    root_nums = []
    for i in range(len(roots)):
        root_nums.append(m('x/.roots[[' + str(i+1) + ']]')) # get the roots in a list, just of the values
        
    root_nums = [abs(root.sage()) for root in root_nums] # change them to sage objects and take their absolute values
    return max(root_nums) # and return the largest


# In[32]:


def get_minimal_direction(sig, method = 'fox'): # doesn't work for dimension 1 homology, but doesn't need to because the minimal direction is always just along the x-axis due to 1-dimensionality.
    '''Returns a list of vectors which are the
    candidates for the vector in the fibered cone
    pointing in the direction of the minimal
    normalized dilatation. COULD IT EVER RETURN MORE THAN 1??
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        tuple
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    tri, angle = taut.isosig_to_tri_angle(sig)
    variables = get_hom_vars(sig)
    f = str(get_Thurston_norm(sig))
    g = str(get_spec_at_e(sig, method))
    g = g.replace('e', 'E') # for mathematica
    var_string = ''
    for variable in variables: # need to turn the variables into a string for mathematica
        if variables.index(variable) != len(variables) - 1:
            var_string += str(variable) + ', '
        else:
            var_string += str(variable)
    
    # Now we solve a Lagrange multipliers problem: minimze the Thurston norm subject to lying on the level set where the normalized dilatation is just equal to the Thurston norm.
    sols = m('sols = NSolve[{Grad[' + f + ', {' + var_string + '}] == L*Grad[' + g + ', {' + var_string + '}], ' + g + ' == 0}, {' + var_string + ',L}, Reals]')
    
    out = []
    for i in range(len(sols)):
        test_point = [m('a' + str(j) + '/.sols[[' + str(i+1) + ']]') for j in range(len(variables))]
        if get_fibered_cone(sig).interior_contains(tuple(test_point)):
            out.append(tuple(test_point))
    return out


# In[33]:


get_minimal_direction(sig)


# In[30]:


def get_norm_dila_log_mathematica(sig, point, method = 'fox'):
    '''Returns the normalized dilatation (in
    the log form, not the exponential form)
    of a homology class.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    point: tuple
        The point in H_2(M, dM; Z) at which you want to
        get the normalized dilatation.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    
    Return type:
        sage.symbolic.expression.Expression
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    keys = get_hom_vars(sig)
    
    if len(keys) != len(point):
        return 'This point is not of the correct dimension.'
    
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitute = {} # Need to make a dictionary to keep track of which variables we are substituting for which variable, since the get_spec function returns the specialization with variables and we need to plug in our point.
    for key, coord in zip(keys, point):
        substitute[key] = coord
    
    return get_Thurston_norm(sig).subs(substitute)*math.log(get_dila_mathematica(sig, point, method))


# In[31]:


type(get_norm_dila_log_mathematica(sig, (-1, 2.9, 1)))


# In[32]:


def get_min_norm_dila_log_mathematica_approx(sig, method = 'fox', rounding = 1):
    '''Returns (approximately) the minimal normalized
    dilatation for the given fibered face of the given 3-manifold.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    method: str, optional
        The method (from the Veering repository taut_polynomial file)
        that we use to compute the taut_polynomial.
        Default: 'fox'
        Options: 'original', 'fox', 'tree' (fox and tree tend to be fastest I think)
    rounding: int, optional
        The number of decimal places to which you want to round
        the entries of the minimal direction to speed up computation.
        Default: 1.
    
    Return type:
        sage.symbolic.expression.Expression
    '''
    if method not in ['original', 'fox', 'tree']:
        raise Exception('Not a valid choice for method.')
        
    tri, angle = taut.isosig_to_tri_angle(sig)
    if tri.homology(1).rank() == 1:
        if get_fibered_cone(sig).contains(1,):
            return get_norm_dila_log_mathematica(sig, (1,), method) # These fibered cones always point along the positive or negative x-axis.
        else:
            return get_norm_dila_log_mathematica(sig, (-1,), method)
    min_direction = get_minimal_direction(sig, method)[0] # REMOVE THE [0] IF GET_MINIMAL_DIRECTION CHANGES OUTPUT TO JUST ONE POINT!
    point = tuple(round(direction, rounding) for direction in min_direction) # round for speed purposes
    return get_norm_dila_log_mathematica(sig, point, method)


# In[33]:


get_min_norm_dila_log_mathematica_approx(sig)


# In[34]:


def get_boundary_triangulations(sig):
    '''Returns a list of ordered pairs.
    Each element of the list corresponds to a boundary
    component of the 3-manifold. The first entry in the
    ordered pair gives a Regina triangulation for that
    component, and the second entry gives how it is
    labeled. (see Regina's buildLink() and
    buildLinkInclusion() for details).
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[tuple(regina.engine.Triangulation2, regina.engine.Isomorphism3)]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    return [(vertex.buildLink(),vertex.buildLinkInclusion()) for vertex in tri.vertices()]


# In[35]:


get_boundary_triangulations(sig)


# In[36]:


def get_C2dM_to_C1dM(sig):
    '''Returns a list of matrices; one for each
    boundary component of the triangulation giving
    the boundary map from 2-chains to 1-chains.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[sage.matrix.matrix_integer_dense.Matrix_integer_dense]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    retval = []
    for vertex in tri.vertices():
        bdry_comp = vertex.buildLink()
        bdry_map = bdry_comp.boundaryMap(2)
        retval.append(matrix([[bdry_map.entry(j,i).pythonValue() for i in range(bdry_map.columns())] for j in range(bdry_map.rows())]))
    
    return retval


# In[37]:


get_C2dM_to_C1dM(sig)


# In[38]:


def get_C1dM_to_C0dM(sig):
    '''Returns a list of matrices; one for each
    boundary component of the triangulation giving
    the boundary map from 1-chains to 0-chains.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[sage.matrix.matrix_integer_dense.Matrix_integer_dense]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    retval = []
    for vertex in tri.vertices():
        bdry_comp = vertex.buildLink()
        bdry_map = bdry_comp.boundaryMap(1)
        retval.append(matrix([[bdry_map.entry(j,i).pythonValue() for i in range(bdry_map.columns())] for j in range(bdry_map.rows())]))
    
    return retval


# In[39]:


get_C1dM_to_C0dM(sig)


# In[40]:


def get_boundary_homology_with_gens(sig):
    '''Returns a list of homology groups with
    generators; one for each boundary component.
    Having the generators will help determine the
    boundary map from the LES in homology.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[list[tuple(sage.homology.homology_group.HomologyGroup_class_with_category, sage.homology.chain_complex.ChainComplex_class_with_category.element_class)]]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    retval = []
    for i in range(tri.countBoundaryComponents()):
        chain_complex = sage.homology.chain_complex.ChainComplex({2: get_C2dM_to_C1dM(sig)[i], 1: get_C1dM_to_C0dM(sig)[i]}, degree_of_differential=-1)
        retval.append(chain_complex.homology(1, generators = True))
                      
    return retval


# In[41]:


get_boundary_homology_with_gens(sig)


# In[42]:


def get_prongs_matrix(sig):
    '''Returns a matrix which, when multiplying a
    2-chain by this matrix, gives the total number
    of prongs on each of the boundary components.
    To get the number of prongs on any given
    singularity, note that all singularities
    coming from a given boundary component have the
    same number of prongs. So we simply divide the
    number of singularities on a given boundary
    component by the total number of prongs on that
    boundary component. This function works by looking
    at Anna Parlak's lower track. Specifically for each
    triangle, look at its bottom embedding. The lower
    track contributes a prong on the vertex shared with
    the top edge (see pictures in her paper). So we loop
    over all triangles and find the vertex to which the
    lower track contributes a prong.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        sage.matrix.matrix_integer_dense.Matrix_integer_dense
    ''' 
    tri, angle = taut.isosig_to_tri_angle(sig)
    
    mat = [[0 for i in range(tri.countVertices())] for j in range(tri.countTriangles())] # initialize our matrix

    bottom_embeddings = transverse_taut.top_bottom_embeddings_of_faces(tri, angle)[1] # the top embeddings of all faces
    for triangle in tri.triangles(): # loop over each face
        bottom_embed = bottom_embeddings[triangle.index()] # and find the bottom embedding of that face
        tet = bottom_embed.tetrahedron() # note the tetrahedron for the bottom embedding
    
        tet_vert_coors = transverse_taut.is_transverse_taut(tri, angle, return_type = "tet_vert_coorientations") # need this to find the top vertices of the given tet
        tet_top_verts = transverse_taut.get_tet_top_vert_nums(tet_vert_coors, tet.index()) # find the vertices of the top edge of the given tet
    
        for i in range(3): # now find the vertex on the top edge which is also a vertex on our bottom face
            current_vert = bottom_embed.vertices()[i]
            if current_vert in tet_top_verts:
                tet_vert = current_vert
            
        vert = tet.vertex(tet_vert).index() # and get its index in the triangulation of the whole 3-manifold
        mat[triangle.index()][vert] = 1 # and add one to the row corresponding to this triangle and the column corresponding to the given vertex
        
    return matrix(mat).transpose()


# In[43]:


get_prongs_matrix(sig)


# In[44]:


def get_boundary_tri_labelings_from_tet_and_vert(sig):
    '''This is effectively an auxiliary function
    for get_C2M_to_C1dM. Returns a dictionary where
    keys are tuples giving a tetrahedron number and
    a vertex on that tetrahedron, and associated value
    is a tuple where the first entry is the boundary
    component number on which that vertex lies and the
    second entry is the index of the associated triangle
    in the triangulation of that boundary component.
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        dict{tuple(int, int) : tuple(int, int)}
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    assocs = {}
    bndries = [(vertex.buildLink(), vertex.buildLinkInclusion()) for vertex in tri.vertices()]
    for bndry in bndries:
        bndry_comp = bndries.index(bndry)
        for triangle in bndry[0].triangles():
            triangle_ind = triangle.index()
            tet = bndry[1].tetImage(triangle_ind) # get the tetrahedron it comes from
            vert = bndry[1].facePerm(triangle_ind)[3] # get the vertex on the tetrahedron that it comes from
            assocs[(tet,vert)] = (bndry_comp, triangle_ind)
            
    return assocs


# In[45]:


get_boundary_tri_labelings_from_tet_and_vert(sig)


# In[46]:


def get_C2M_to_C1dM(sig):
    '''Returns a list of matrices. Each matrix
    corresponds to one boundary component. They
    each give the map assigning triangles in the
    triangulation of the 3-manifold to the edges
    in each respective boundary component coming
    from that triangle. The signs are in agreement
    with the veering orientations. To see how this
    function was constructed, it helps to look at
    pictures in Dilas Week 2 on notability.<<<<<< UPLOAD PICTURES
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[sage.matrix.matrix_integer_dense.Matrix_integer_dense]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    bndries = [(vertex.buildLink(),vertex.buildLinkInclusion()) for vertex in tri.vertices()] # build each boundary component triangulation and the data encoding its inclusion
    assocs = get_boundary_tri_labelings_from_tet_and_vert(sig) # get the dictionary with keys = (tet, vertex) and associated value (boundary component, triangle number)
    #print(assocs) # For testing
    maps = [[[0 for i in range(bndries[k][0].countEdges())] for j in range(tri.countTriangles())] for k in range(tri.countVertices())] # initialize the appropriate number of appropriate dimension matrices
    #print(maps) # For testing
    
    for l in range(tri.countTriangles()): # loop over all the triangles in the 3-manifold
        #triangle = tri.triangles()[l] # store the triangle itself
        top_embed = transverse_taut.top_bottom_embeddings_of_faces(tri, angle)[0][l] # we will work with its top embedding
        tet = top_embed.tetrahedron() # the tetrahedron it belongs to
        tet_ind = tet.index()
        
        alpha = top_embed.vertices()[0] # this is the zeroth vertex of the triangle (as labeled by a vertex on tet)
        beta = top_embed.vertices()[1] # this is the first
        gamma = top_embed.vertices()[2] # this is the second
        
        (bndry_alpha,triangle_alpha) = assocs[(tet_ind,alpha)] # this is the boundary component and the triangle number that corresponds to the vertex alpha
        (bndry_beta,triangle_beta) = assocs[(tet_ind,beta)] # and so on
        (bndry_gamma,triangle_gamma) = assocs[(tet_ind,gamma)]
        #print(bndry_alpha, bndry_beta,bndry_gamma)
        
        x_alpha = bndries[bndry_alpha][1].facePerm(triangle_alpha).inverse()[gamma] # this is the first corner of the alpha triangle (as labeled in the triangulation of the boundary) when traversing according to the RHR alpha-->beta-->gamma
        y_alpha = bndries[bndry_alpha][1].facePerm(triangle_alpha).inverse()[beta] # this is the second
        x_beta = bndries[bndry_beta][1].facePerm(triangle_beta).inverse()[alpha] # and so on
        y_beta = bndries[bndry_beta][1].facePerm(triangle_beta).inverse()[gamma]
        x_gamma = bndries[bndry_gamma][1].facePerm(triangle_gamma).inverse()[beta]
        y_gamma = bndries[bndry_gamma][1].facePerm(triangle_gamma).inverse()[alpha]
        
        #print(str(i) + ' (' + str(x_alpha) + str(y_alpha) + ')') # For testing
        for edge in bndries[bndry_alpha][0].edges(): # loop over the edges in the boundary component containing the alpha vertex
            ind = edge.index()
            for embed in edge.embeddings(): # loop over the embeddings
                if str(embed) == str(triangle_alpha) + ' (' + str(x_alpha) + str(y_alpha) + ')': # if one embedding matches the vertex ordering x_alpha --> y_alpha (within triangle i)
                    maps[bndry_alpha][l][ind] = 1 # then we traversed this edge in the positive direction
                elif str(embed) == str(triangle_alpha) + ' (' + str(y_alpha) + str(x_alpha) + ')': # if one embedding matches the vertex ordering y_alpha --> x_alpha (within triangle i)
                    maps[bndry_alpha][l][ind] = -1 # then we traversed this edge in the negative direction
                    
        # And repeat for the other edges:
        for edge in bndries[bndry_beta][0].edges():
            ind = edge.index()
            for embed in edge.embeddings():
                if str(embed) == str(triangle_beta) + ' (' + str(x_beta) + str(y_beta) + ')':
                    maps[bndry_beta][l][ind] = 1
                elif str(embed) == str(triangle_beta) + ' (' + str(y_beta) + str(x_beta) + ')':
                    maps[bndry_beta][l][ind] = -1
                    
        for edge in bndries[bndry_gamma][0].edges():
            ind = edge.index()
            for embed in edge.embeddings():
                if str(embed) == str(triangle_gamma) + ' (' + str(x_gamma) + str(y_gamma) + ')':
                    maps[bndry_gamma][l][ind] = 1
                elif str(embed) == str(triangle_gamma) + ' (' + str(y_gamma) + str(x_gamma) + ')':
                    maps[bndry_gamma][l][ind] = -1
                    
    # As of now, the traversals are happening according to the RHR wrt to the labelings.
    # We now fix things to they agree with the veering coorientation.
    #print(maps)
    agreement = transverse_taut.is_transverse_taut(tri, angle, return_type = "face_coorientations") # This gives a list telling us whether the RHR orientations of the regina numbering agrees or disagrees with the veering coorientation. (See Veering repo for details)
    for p in range(len(bndries)):
        for q in range(tri.countTriangles()):
            #print(maps[p][q])
            maps[p][q] = [maps[p][q][r]*agreement[q] for r in range(len(maps[p][q]))]
            #print(maps[p][q])
                
    retval = [matrix(maps[r]).transpose() for r in range(len(bndries))]
        
    return retval


# In[47]:


get_C2M_to_C1dM(sig)


# In[48]:


def get_LES_boundary_maps(sig):
    '''Returns a list (one entry for each
    boundary component) of doubles. Each
    double represents the image of an
    arbitrary homology class under the
    boundary map for the LES sequence in
    homology H_2(M) --> H_1(dM).
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[tuple(sage.symbolic.expression.Expression)]
    '''
    tri, angle = taut.isosig_to_tri_angle(sig)
    retval = []
    for i in range(tri.countVertices()): # loop over the boundary components
        gens = [get_boundary_homology_with_gens(sig)[i][j][1].vector(1) for j in range(2)] # these are the 1-chains whose images in homology generate the first homology of this boundary component.
        #print(gens)
        
        # What we want to do is write an arbitrary vector (vector with variables) in C_1(dM) which represents
        # a homology class (that is, it is in the kernel of the map C_1(dM) --> C_0(dM)) as a linear combination
        # of the generators above and elements in the image of C_2(dM) --> C_1(dM), so that the class in homology
        # will simply be the coefficients of the generator vector.
        
        C2dM_to_C1dM_comp_i = get_C2dM_to_C1dM(sig)[i] # get boundary map from 2-chains to 1-chains in this boundary component so that we can find the vectors which generate the image.
        basis_for_image_of_boundary21 = [C2dM_to_C1dM_comp_i.column(j) for j in C2dM_to_C1dM_comp_i.pivots()] # this is a basis for the image
        #print(basis_for_image_of_boundary21)
        preimage = get_2_chain_from_point_in_homology(sig) # given an arbitrary point in homology, this is a 2-chain representing that homology class.
        C2M_to_C1dM_comp_i = get_C2M_to_C1dM(sig)[i] # this is the boundary map (on the chain level) that we are trying to get on the homology level
        arb_image_in_C1dM_comp_i = C2M_to_C1dM_comp_i*preimage # this is the vector we want to write as a linear combo of the basis and gens above
        
        mat0 = [basis_for_image_of_boundary21[i] for i in range(len(basis_for_image_of_boundary21))]
        mat0.insert(0, gens[1])
        mat0.insert(0, gens[0])
        mat0 = matrix(mat0).transpose()
        #print(mat0)
        sols = mat0.solve_right(arb_image_in_C1dM_comp_i)
        sols0 = (sols[0],sols[1])
        #print(matrix(sols0))
        retval.append(sols0)
    return retval


# In[49]:


get_LES_boundary_maps(sig)


# In[50]:


def get_LES_boundary_matrices(sig):
    '''Returns a list (one entry for each
    boundary component) of matrices. Each
    matrix represents the the boundary map
    for the LES sequence in homology
    H_2(M) --> H_1(dM).
    
    Parameters
    __________
    sig: str
        The isomorphism signature (together with 
        veering structure) of a veering triangulation.
        That is, an entry in the veering census.
    
    Return type:
        list[sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense]
    ''' 
    tri, angle = taut.isosig_to_tri_angle(sig)
    retval = []
    for i in range(tri.countVertices()): # loop over the boundary components
        gens = [get_boundary_homology_with_gens(sig)[i][j][1].vector(1) for j in range(2)] # these are the 1-chains whose images in homology generate the first homology of this boundary component.
        #print(gens)
        
        # What we want to do is write an arbitrary vector (vector with variables) in C_1(dM) which represents
        # a homology class (that is, it is in the kernel of the map C_1(dM) --> C_0(dM)) as a linear combination
        # of the generators above and elements in the image of C_2(dM) --> C_1(dM), so that the class in homology
        # will simply be the coefficients of the generator vector.
        
        C2dM_to_C1dM_comp_i = get_C2dM_to_C1dM(sig)[i] # get boundary map from 2-chains to 1-chains in this boundary component so that we can find the vectors which generate the image.
        basis_for_image_of_boundary21 = [C2dM_to_C1dM_comp_i.column(j) for j in C2dM_to_C1dM_comp_i.pivots()] # this is a basis for the image
        #print(basis_for_image_of_boundary21)
        preimage = get_2_chain_from_point_in_homology(sig) # given an arbitrary point in homology, this is a 2-chain representing that homology class.
        C2M_to_C1dM_comp_i = get_C2M_to_C1dM(sig)[i] # this is the boundary map (on the chain level) that we are trying to get on the homology level
        arb_image_in_C1dM_comp_i = C2M_to_C1dM_comp_i*preimage # this is the vector we want to write as a linear combo of the basis and gens above
        
        mat0 = [basis_for_image_of_boundary21[i] for i in range(len(basis_for_image_of_boundary21))]
        mat0.insert(0, gens[1])
        mat0.insert(0, gens[0])
        mat0 = matrix(mat0).transpose()
        #print(mat0)
        sols = mat0.solve_right(arb_image_in_C1dM_comp_i)
        sols0 = (sols[0],sols[1])
        #print(matrix(sols0))
        retmat = [[0 for j in range(tri.homology().rank())] for k in range(2)]
        for j in range(2):
            for k in range(len(retmat[0])):
                retmat[j][k] = sols0[j].coefficient(var('a' + str(k)))
        retval.append(matrix(retmat))
    return retval


# In[51]:


get_LES_boundary_matrices(sig)


# In[52]:


#'''This was just to test that my get_C2M_to_C1dM function was working correctly. Used the KKT example.'''
#
#testComp = tri.vertices()[0].buildLink()
#testCompData = tri.vertices()[0].buildLinkInclusion()
#
#for i in range(testComp.countTriangles()):
#    print('Triangle ' + str(i) + ' lies on the ' + str(testCompData.tetImage(i)) + ' tetrahedron')
#    print('on the ' + str(testCompData.facePerm(i)[3]) + ' vertex')
#    print('with its 0 vertex connected to ' + str(testCompData.facePerm(i)[0]))
#    print('its 1 vertex connected to ' + str(testCompData.facePerm(i)[1]))
#    print('and its 2 vertex connected to ' + str(testCompData.facePerm(i)[2]) + ' of the tetrahedron')


# # For working with specific examples:

# Note that this is not the fastest way to deal with specific examples since every function recalculates everything it needs to.

# In[ ]:


def get_num_punctures(sig, point):
    '''This function returns a list. The first entry
    of the list is another list, where each entry
    corresponds to a boundary component. The numbers
    give the number of singularities of the fiber
    corresponding to the given point which lie on the
    corresponding boundary component. The last number
    (the second entry of the outer list) gives the total
    number of singularities.'''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    assert len(point) == tri.homology().rank(), "The dimension of the point you entered and the second homology are not the same rank."
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitution_dict = {get_hom_vars(sig)[i]:point[i] for i in range(len(point))} # substitute the given point for our variable place holders.
    #print(substitution_dict)
    boundary_maps = get_LES_boundary_maps(sig)
    retval = []
    for comp in boundary_maps:
        retval.append(gcd(comp[0].subs(substitution_dict),comp[1].subs(substitution_dict))) # the number of singularities is the gcd of the images under the LES boundary maps
    
    return [retval, sum(retval)]


# In[ ]:


def get_genus(sig, point):
    '''This function returns the genus
    of the surface fiber associated to the
    given point in the homology of the given
    3-manifold.'''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    assert len(point) == tri.homology().rank(), "The dimension of the point you entered and the second homology are not the same rank."
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("Either the point you entered is not in the interior of the fibered cone, or this is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitution_dict = {get_hom_vars(sig)[i]:point[i] for i in range(len(point))} # substitute the given point for our variable place holders.
    t_norm = get_Thurston_norm(sig).subs(substitution_dict)
    num_punctures = get_num_punctures(sig, point)[1]
    
    # Just follow the formula for the Thurston norm.
    return (t_norm + 2 - num_punctures)/2


# In[ ]:


def get_num_prongs(sig, point): # DON'T THINK THIS WORKS WITH VARIABLES
    '''This function returns a list of lists.
    One list for each boundary component. The
    contents of the list is the number of prongs
    on a singularity lying on that boundary component,
    repeated the appropriate number of times so as to
    agree with the number of singularities on this
    boundary component. Of course, these numbers are
    describing the surface fiber associated to the
    given point in the homology of the given
    3-manifold.'''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    assert len(point) == tri.homology().rank(), "The dimension of the point you entered and the second homology are not the same rank."
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("This is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    substitution_dict = {get_hom_vars(sig)[i]:point[i] for i in range(len(point))} # substitute the given point for our variable place holders.
    two_chain_with_vars = get_2_chain_from_point_in_homology(sig).subs(substitution_dict) # get the corresponding 2-chain and substitute our point for the relevant variables
    no_vars = [] # two_chain_with_vars still contains the free variables. We want to set those all to 0.
    for term in two_chain_with_vars:
        no_vars.append(term.subs({v:0 for v in term.variables()}))

    no_vars = vector(no_vars)
    total_prongs_each_boundary_comp = get_prongs_matrix(sig)*no_vars # multiply our 2-chain by the prongs matrix defined above. This gives the total number of prongs on each boundary component.
    punctures = get_num_punctures(sig, point)
    return [[total_prongs_each_boundary_comp[i]/punctures[0][i]]*int(punctures[0][i]) for i in range(tri.countBoundaryComponents())]


# In[ ]:


# Taken from:
# https://ask.sagemath.org/question/62717/fastest-way-to-solve-non-negative-linear-diophantine-equations/

from sage.numerical.mip import MIPSolverException

def mysolve(A,V):
    '''This is an auxiliary function used for
    finding only nonnegative solutions to linear
    equations which we need for testing the functions
    above. We use it to find nonnegative carried 2-chains
    so we can use the veering functions to build the
    surface fiber and compare the relvant data.'''
    
    n = len(A)      # number of equations
    assert n == len(V)
    m = len(A[0])   # number of variables

    milp = MixedIntegerLinearProgram(solver='ppl')
    milp.set_objective( None )
    x = milp.new_variable(integer=True, nonnegative=True)
    for i in range(n):
        milp.add_constraint( sum(A[i][j] * x[j] for j in range(m)) == V[i] )

    try:
        milp.solve()
    except MIPSolverException:
        return []       # no solutions

    X = milp.get_values(x)
    
    return vector([X[i] for i in range(m)])


# In[ ]:


def get_carried_2_chain(sig, point):
    '''This is an auxiliary function used
    to find nonnegative carried 2-chains so
    we can use the veering functions to build
    the surface fiber and compare the relvant
    data.'''
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    assert len(point) == tri.homology().rank(), "The dimension of the point you entered and the second homology are not the same rank."
    try:
        cone = get_fibered_cone(sig)
        assert cone.interior_contains(point), "The point is not in the interior of the fibered cone."
    except:
        print("This is an example where the fibered cone is not built into the Veering taut_polytope file (not sure why). Be carefule to make sure your point is in the interior of the fibered cone.")
    
    P = taut_polytope.projection_to_homology(sig)
    Pmatrix = list(P)
    Pmatrix = [list(row) for row in Pmatrix] # the auxiliary function takes in matrices as lists of lists
    C2M_to_C1Mmatrix = list(get_C2M_to_C1M(sig))
    C2M_to_C1Mmatrix = [list(row) for row in C2M_to_C1Mmatrix]
    A = Pmatrix + C2M_to_C1Mmatrix # We want our vector to be both in the preimage of (a,b,c) and in the kernel of C2M_to_C1M so we just stack the matrices.
    V = vector([point[i] for i in range(len(point))] + [0 for j in range(tri.countEdges())])
    #V = vector([point[0],point[1],point[2],0,0,0,0,0,0,0,0]) # And modify the target vector accordingly.

    return mysolve(A,V)


# In[ ]:


def test_a_point(sig, point):
    '''For testing.'''
    tri, angle = taut.isosig_to_tri_angle(sig)
    print('calculated genus: ', get_genus(sig, point))
    print('calculated punctures: ', get_num_punctures(sig, point))
    print('calculated prongs: ', get_num_prongs(sig, point))
    print('veering data ((genus, punctures), [prongs]): ', carried_surface.stratum(tri, angle, get_carried_2_chain(sig, point)))


# In[ ]:




