#!/usr/bin/env python
# coding: utf-8

# In[2]:


from veering import taut_polytope
from veering import taut
from veering import taut_homology
from veering import taut_carried
from veering import fundamental_domain
from veering import taut_polynomial
from veering import transverse_taut
from sage.plot.contour_plot import ContourPlot
from itertools import combinations
from itertools import permutations
from veering import carried_surface
import time

m = mathematica


# In[3]:


def get_extreme_surfs(sig):
    # This function gets the branch equation form of the surfaces which span the fibered cone.
    
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


# In[4]:


def get_corners_of_fibered_face(sig):
    # This function gets the corners of the fibered face determined by the given veering triangulation.
    
    extreme_rays = taut_polytope.cone_in_homology(sig) # These are rays spanning the cone of carried classes
    extreme_surfs = get_extreme_surfs(sig)
    sums = [] # we add up the number of triangles in the boundary surfaces
    for surf in extreme_surfs:
        count = 0
        for i in range(len(surf)):
            count = count + surf[i]
        sums.append(count)

    tnorms = [term/2 for term in sums] # their Thurston norms are half the number of boundary components
        
    corners = []
    for i in range(len(extreme_rays)):
        corners.append([0]*len(extreme_rays[i]))
    
    for i in range(len(corners)): # the corners of the fibered face are the boundaries divided by the Thurston norm.
        for j in range(len(corners[i])):
            corners[i][j] = extreme_rays[i][j]/tnorms[i]
        
    return corners


# In[5]:


def get_fibered_cone(sig):
    return Cone(get_corners_of_fibered_face(sig))


# In[6]:


def get_fibered_cone_plot(sig):
    return get_fibered_cone(sig).plot(mode='box', xmin = -10, xmax = 10, ymin = -10, ymax = 10, zmin = -10, zmax = 10, ray_thickness=2, point_size=10)


# In[7]:


def get_surf_from_point_in_homology(sig, point): #CHECK
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    if tri.homology().rank() == 1:
        return point[0]*get_extreme_surfs(sig)[0]
    
    extreme_surfs = get_extreme_surfs(sig)
    P = taut_polytope.projection_to_homology(sig) # This is the matrix to multiply by to get branched surface tuple to a second homology class.
    
    variables = [var('x' + str(i)) for i in range(len(extreme_surfs))]
    
    eqn_pieces = [variables[i]*P*extreme_surfs[i] for i in range(len(extreme_surfs))]
        
    term = 0
    for k in range(len(eqn_pieces)):
        term += eqn_pieces[k]
    
    eqns = [term[i] == point[i] for i in range(len(point))]
    
    lin_comb = solve(eqns, variables)
    
    retval = 0
    for i in range(len(lin_comb[0])):
        retval += lin_comb[0][i].rhs()*extreme_surfs[i]

    return retval

#get_surf_from_point_in_homology('fLLQcbecdeepuwsua_20102', (var('a'), var('b')))


# In[8]:


def Thurston_norm_nums(sig, point):
    extreme_rays = taut_polytope.cone_in_homology(sig) # These are rays spanning the cone of carried classes
    
    if not get_fibered_cone(sig).interior_contains(point):
        raise ValueError('Your point is not in the fibered cone specified by the given veering triangulation')
    
    if len(point) != len(extreme_rays[0]):
        raise ValueError('Your point must be of dimension', len(extreme_rays[0])) 

    surf = get_surf_from_point_in_homology(sig, point)
    prenorm = 0
    for tri in surf:
        prenorm = prenorm + tri
        
    return abs(prenorm/2) #abs b/c for some reason some Betti number 1 cases were returning negative - the ones whose cones point along the negative x-axis.

#Thurston_norm_nums('fLLQcbecdeepuwsua_20102',(-1,1))


# In[9]:


def Thurston_norm_vars(sig, point): # This version just doesn't check if the input is in the cone.
    extreme_rays = taut_polytope.cone_in_homology(sig) # These are rays spanning the cone of carried classes
    
    if len(point) != len(extreme_rays[0]):
        raise ValueError('Your point must be of dimension', len(extreme_rays[0])) 
        
    surf = get_surf_from_point_in_homology(sig, point)
    prenorm = 0
    for tri in surf:
        prenorm = prenorm + tri
        
    return prenorm/2


# In[10]:


def get_spec(sig, point):
    # returns the specialization of the taut polynomial at a point in second homology.
    
    tri, angle = taut.isosig_to_tri_angle(sig)
    poly = taut_polynomial.taut_polynomial(tri, angle)
    unspec_monos = poly.monomials()
    monos = [] # the monomials in the specialization
    for i in range(poly.number_of_terms()):
        monos.append(unspec_monos[i].lc()*x^(sum([poly.exponents()[i][j]*point[j] for j in range(len(point))])))
        
    spec = 0
    for i in range(len(monos)):
        spec = spec + monos[i]*poly.monomial_coefficient(poly.monomials()[i])
        
    return spec

#get_spec('gvLQQcdeffeffffaafa_201102', (var('a'), var('b'), var('c')))


# In[11]:


def get_dila(sig, point):
    spec = get_spec(sig, point)    
    sols = solve(spec == 0, x, to_poly_solve=True)
    
    sols2 = []
    for i in range(len(sols)):
        sols2.append(sols[i].rhs())
       
    for j in range(len(sols2)):
        sols2[j] = abs(sols2[j])
    
    return max(sols2)

#get_dila('gvLQQcdeffeffffaafa_201102', (2,6,1))


# In[12]:


def get_norm_dila_log(sig, point):
    return Thurston_norm_nums(sig, point)*math.log(get_dila(sig, point))

#get_norm_dila_log('gvLQQcdeffeffffaafa_201102', (2,6,1))


# In[13]:


def get_spec_at_e(sig, point): # No need to use for Betti number 1 examples
    return get_spec(sig, point)(x = e)

#get_spec_at_e('ivLLQQccdgffhghhvvaaaaavv_01111220', (var('a'), var('b'), var('c')))


# In[14]:


def get_fibered_cone_and_levset_plot(sig):
    a,b,c = var('a,b,c')
    levset = implicit_plot3d(get_spec_at_e(sig, (a,b,c)) == 0, (-10,10), (-10,10), (-10,10), cmap=['green'])
    Plot = levset + get_fibered_cone_plot(sig)
    return Plot


# In[15]:


def get_minimal_direction(sig): # No need to use for Betti number 1 examples
    tri, angle = taut.isosig_to_tri_angle(sig)
    variables = [var('a' + str(i)) for i in range(tri.homology().rank())]
    f = str(Thurston_norm_vars(sig, variables))
    g = str(get_spec_at_e(sig, variables))
    g = g.replace('e', 'E')
    var_string = ''
    for variable in variables:
        if variables.index(variable) != len(variables) - 1:
            var_string += str(variable) + ', '
        else:
            var_string += str(variable)
            
    sols = m('sols = NSolve[{Grad[' + f + ', {' + var_string + '}] == L*Grad[' + g + ', {' + var_string + '}], ' + g + ' == 0}, {' + var_string + ',L}, Reals]')
    
    for i in range(len(sols)):
        test_point = [m('a' + str(j) + '/.sols[[' + str(i+1) + ']]') for j in range(len(variables))]
        if get_fibered_cone(sig).interior_contains(tuple(test_point)):
            return tuple(test_point)

#get_minimal_direction('gvLQQcdeffeffffaafa_201102')


# In[16]:


def get_dila_mathematica(sig, point):
    poly = str(get_spec(sig, point))
    poly = poly.replace("sqrt(x)", "Sqrt[x]")
    #print(poly)
    roots = m('roots = NSolve[' + poly + ' == 0, x]')
    #print(roots)
    root_nums = []
    for i in range(len(roots)):
        root_nums.append(m('x/.roots[[' + str(i+1) + ']]'))
        
    root_nums = [abs(root.sage()) for root in root_nums]
    return max(root_nums)


# In[17]:


def get_norm_dila_log_mathematica(sig, point):
    return Thurston_norm_nums(sig, point)*math.log(get_dila_mathematica(sig, point))


# In[18]:


def get_min_norm_dila_log_mathematica_approx(sig):
    tri, angle = taut.isosig_to_tri_angle(sig)
    if tri.homology().rank() == 1:
        if get_fibered_cone(sig).contains(1,):
            return get_norm_dila_log_mathematica(sig, (1,)) # These fibered cones always point along the positive or negative x-axis.
        else:
            return get_norm_dila_log_mathematica(sig, (-1,))
    min_direction = get_minimal_direction(sig)
    point = tuple(round(direction, 1) for direction in min_direction)
    return get_norm_dila_log_mathematica(sig, point)


# In[19]:


'''with open('MinimalDilatationExamplesWithBetti1.txt', 'w') as g:
    with open('veering_census.txt') as f:
        lines = f.readlines()

    for i in range(len(lines)):
        lines[i] = lines[i][:-1]
    
    for sig in lines:
        tri, angle = taut.isosig_to_tri_angle(sig)
        try:
            start = time.time()
            min_norm_dila = str(get_min_norm_dila_log_mathematica_approx(sig))
            end = time.time()
            g.write('sig: ' + sig + '\t rank: ' + str(tri.homology().rank()) + '\t min_norm_dila: ' + min_norm_dila + '\t time: ' + str(end-start))
            g.write('\n')
        except Exception: 
            g.write(sig + ' had an error. \t rank: ' + str(tri.homology().rank()))
            g.write('\n')
            pass'''


# In[ ]:




