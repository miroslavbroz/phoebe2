/* 
  Wrappers for calculations concerning the generalized Kopal potential or
  Roche lobes.
   
  
  Need to install for Python.h header:
    apt-get install python-dev
  
  Ref:
  
  Python C-api: 
  * https://docs.python.org/2/c-api/index.html 

  Numpi C-api manual:
  * http://docs.scipy.org/doc/numpy/reference/c-api.html
  * http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
  * http://scipy.github.io/old-wiki/pages/Cookbook/C_Extensions/NumPy_arrays
  
  Wrapping tutorial:
  * http://intermediate-and-advanced-software-carpentry.readthedocs.io/en/latest/c++-wrapping.html
  
  Author: Martin Horvat, June 2016
*/

#include <iostream>
#include <vector>
#include "gen_roche.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
 
// providing the Python interface -I/usr/include/python2.7/
#include <Python.h>
#include <numpy/arrayobject.h>

/*
  Creating PyArray from std::vector.
  
  Ref: 
  http://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.ndarray.html
*/ 
template <class T>
PyObject  *PyArray_SimpleNewFromVector(
  int np, 
  npy_intp *dims, 
  int typenum, 
  void *V){
  
  int size = 1;
  for (int i = 0; i < np; ++i) size *= dims[i];
  size *= sizeof(T);
  
  // Note: p is interpreted in C-style contiguous fashion
  void *p = malloc(size);
  memcpy(p, V, size);
  
  return PyArray_SimpleNewFromData(np, dims, typenum, p); 
}


/*
  Python wrapper for C++ code:
  
  Calculate critical potential, that is the value of the generalized 
  Kopal potential in Lagrange points.
  
  Python:
    
    omega = roche_critical_potential(q, F, d)
  
  where parameters are floats
  
    q = M2/M1 - mass ratio
    F - synchronicity parameter
    d - separation between the two objects
  
  and returns
  
    omega - tuple of 3 float numbers
*/

static PyObject *roche_critical_potential(PyObject *self, PyObject *args) {
    
    // parse input arguments   
    double q, F, delta, omega[3];
    
    if (!PyArg_ParseTuple(args, "ddd", &q, &F, &delta))
        return NULL;
        
    // calculate critical potential
    gen_roche::critical_potential(omega, q, F, delta);
    
    // create a tuple and store omega[3]
    PyObject *omega_o = PyTuple_New(3);
    
    for (int i = 0; i < 3; ++i) 
      PyTuple_SetItem(omega_o, i, PyFloat_FromDouble(omega[i]));
      
    return omega_o;
}


/*
  Python wrapper for C++ code:
  
  Calculate height h the left (right) Roche lobe at the position 
    
    (0,0,h)   (  (d,0,h) )
    
  of the primary (secondary) star. Roche lobe(s) is defined as equipotential 
  of the generalized Kopal potential Omega:

      Omega_0 = Omega(x,y,z)
  
  Python:
    
    h = roche_Rpole(q, F, d, Omega0, choice)
  
  where parameters are floats
  
    q = M2/M1 - mass ratio
    F - synchronicity parameter
    d - separation between the two objects
    Omega - value potential 
    choice - 0 for discussing left lobe, 1 for discussing right lobe
  
  and return float
  
    h : height of the lobe's pole
*/


static PyObject *roche_pole(PyObject *self, PyObject *args) {
    
    // parse input arguments   
    int choice;
    
    double q, F, delta, Omega0;
    
    if (!PyArg_ParseTuple(args, "ddddi", &q, &F, &delta, &Omega0, &choice))
        return NULL;
     
    if (choice == 0)  // Left lobe
      return PyFloat_FromDouble(gen_roche::poleL(Omega0, q, F, delta));
    
    // Right lobe  
    return PyFloat_FromDouble(gen_roche::poleR(Omega0, q, F, delta));
}


/*
  Python wrapper for C++ code:

    Marching meshing of Roche lobes implicitely defined by generalized 
    Kopal potential:
    
      Omega_0 = Omega(x,y,z)
    
  Python:

    dict = marching_mesh(q, F, d, Omega0, delta, choice, max_triangles, <keyword>=[true,false], ... )
    
  where parameters are floats
   
      q: float = M2/M1 - mass ratio
      F: float - synchronicity parameter
      d: float - separation between the two objects
      Omega0: float - value of the generalized Kopal potential
      delta: float - size of triangles edges projected to tangent space
      choice: integer type
                0 - primary lobe is exists
                1 - secondary lobe is exists
              for overcontacts choice is 0 or 1
              choice controls where is the begining the triangulation
      max_triangles:integer - maximal number of triangles
        if number of triangles exceeds max_triangles it returns NULL  
  
  Returns 
  
    dictionary
  
  with keywords
  
    vertices:
      V[][3]    - 2-rank numpy array of a pairs of vertices 
    
    vnormals:
      NatV[][3] - 2-rank numpy array of normals at vertices
    
    triangles:
      T[][3]    - 2-rank numpy array of 3 indices of vertices 
                composing triangles of the mesh aka connectivity matrix
    
    tnormals:
      NatT[][3] - 2-rank numpy array of normals of triangles
  
    areas:
      A[]       - 1-rank numpy array of areas of triangles
    
    area:
      area      - area of triangles of mesh
    
    volume:
      volume    - volume of body enclosed by triangular mesh
      
    centers:
      C[][3]    - 2-rank numpy array of central points of triangles
                  central points is  barycentric points projected to 
                  Roche lobes
    cnormals:
      NatC[][3]   - 2-rank numpy array of normals of central points
 

  Typically face-vertex format is (V, T) where
  
  V - vertices
  T - connectivity matrix with indices labeling vertices in 
      counter-clockwise orientation so that normal vector is pointing 
      outward
  
  Refs:
  * for face-vertex format see https://en.wikipedia.org/wiki/Polygon_mesh
  * http://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.ndarray.html
  * http://docs.scipy.org/doc/numpy/reference/c-api.array.html#creating-arrays
  * https://docs.python.org/2.0/ext/parseTupleAndKeywords.html
  * https://docs.python.org/2/c-api/arg.html#c.PyArg_ParseTupleAndKeywords
*/

static PyObject *roche_marching_mesh(PyObject *self, PyObject *args, PyObject *keywds) {
  
  //
  // Reading arguments
  //
 
 char *kwlist[] = {
    (char*)"q",
    (char*)"F",
    (char*)"d",
    (char*)"Omega0",
    (char*)"delta",
    (char*)"choice",
    (char*)"max_triangles",
         
    (char*)"vertices", 
    (char*)"vnormals",
    (char*)"triangles", 
    (char*)"tnormals", 
    (char*)"centers", 
    (char*)"cnormals",
    (char*)"areas",
    (char*)"area",
    (char*)"volume",
    NULL};
  
  double q, F, d, Omega0, delta;   
  
  int choice,
      max_triangles;
      
  bool b_vertices = false, 
       b_vnormals = false, 
       b_triangles = false, 
       b_tnormals = false, 
       b_centers = false,
       b_cnormals = false,
       b_areas = false,
       b_area = false,
       b_volume = false;
      
  if (!PyArg_ParseTupleAndKeywords(
      args, keywds,  "dddddii|iiiiiiiii", kwlist,
      &q, &F, &d, &Omega0, &delta, &choice, &max_triangles,
      &b_vertices, 
      &b_vnormals,
      &b_triangles, 
      &b_tnormals,
      &b_centers,
      &b_cnormals,
      &b_areas,
      &b_area,
      &b_volume 
      ))
    return NULL;
     
  //
  // Storing results in dictioonary
  // https://docs.python.org/2/c-api/dict.html
  
  PyObject *results = PyDict_New();
  
        
  //
  // Points on the x - axis 
  //
  
  std::vector<double> x_points;
    
  gen_roche::points_on_x_axis(x_points, Omega0, q, F, d);

  if (x_points.size() == 0) return NULL;
  
  //
  // Chosing the lobe (left or right or overcontact) 
  //
  
  double x0; 
  
  if (choice == 0) {    // left lobe
    x0 = x_points.front();
    if (x0 > 0) return NULL;
  } else {              // right lobe
    x0 = x_points.back();
    if (x0 < d) return NULL;
  }
  
  //
  //  Marching triangulation of the Roche lobe 
  //
    
  double params[5] = {q, F, d, Omega0, x0};
  
  Tmarching<double, Tgen_roche<double>> march(params);  
  
  std::vector<T3Dpoint<double>> V, NatV;
  std::vector<Ttriangle> Tr; 
  
  if (!march.triangulize(delta, max_triangles, V, NatV, Tr)){
    std::cerr << "There is too much triangles\n";
    return NULL;
  }
  
         
  int Nt = Tr.size(), // number of triangles
      Nv = V.size();  // number of vertices

  npy_intp dims[2];
  
  if (b_vertices) {
    dims[0] = Nv;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "vertices", 
      PyArray_SimpleNewFromVector<double>(2, dims, NPY_DOUBLE, V.data())
    );
  }

  if (b_vnormals) {
    dims[0] = Nv;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "vnormals", 
      PyArray_SimpleNewFromVector<double>(2, dims, NPY_DOUBLE, NatV.data())
    );
  }

  
  if (b_triangles) {
    dims[0] = Nt;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "triangles", 
      PyArray_SimpleNewFromVector<int>(2, dims, NPY_INT, Tr.data())
    );
  }
  
  //
  // Calculte the mesh properties
  //
  int vertex_choice = 0;
  
  double area, volume, 
        *p_area = 0, *p_volume = 0;
  
  std::vector<double> *A = 0; 
  std::vector<T3Dpoint<double>> *NatT = 0;
  
  if (b_areas) A = new std::vector<double>;
  
  if (b_area) p_area = &area;
  
  if (b_tnormals) NatT = new std::vector<T3Dpoint<double>>;
  
  if (b_volume) p_volume = &volume;
   
  mesh_attributes(V, NatV, Tr, A, NatT, p_area, p_volume, vertex_choice, true);

  if (b_areas) {
    dims[0] = Nt;
    PyDict_SetItemString(
      results, 
      "areas",
      PyArray_SimpleNewFromVector<double>(1, dims, NPY_DOUBLE, A->data())
    );
    delete A;  
  }
  
  if (b_area)
    PyDict_SetItemString(
      results, 
      "area",
      PyFloat_FromDouble(area)
    );

  if (b_tnormals) {
    dims[0] = Nt;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "tnormals",
      PyArray_SimpleNewFromVector<double>(2, dims, NPY_DOUBLE, NatT->data())
    );
    delete NatT;  
  }


  if (b_volume)
    PyDict_SetItemString(
      results, 
      "volume",
      PyFloat_FromDouble(volume)
    );


  //
  // Calculte the central points
  // 
  
  std::vector<T3Dpoint<double>> *C = 0, *NatC = 0;
  
  if (b_centers) C = new std::vector<T3Dpoint<double>>;
 
  if (b_cnormals) NatC = new std::vector<T3Dpoint<double>>;
 
  march.central_points(V, Tr, C, NatC);
  
  if (b_centers) {
    dims[0] = Nt;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "centers",
      PyArray_SimpleNewFromVector<double>(2, dims, NPY_DOUBLE, C->data())
    );
    delete C;  
  }


  if (b_cnormals) {
    dims[0] = Nt;
    dims[1] = 3;
    PyDict_SetItemString(
      results, 
      "cnormals",
       PyArray_SimpleNewFromVector<double>(2, dims, NPY_DOUBLE, NatC->data())
    );
    delete NatC;  
  }
  
   
  return results;
}

/*
  Python wrapper for C++ code:

    Calculation of visibility of triangles
    
  Python:

    M = triangle_mesh_visibility(v, V, T, N)
    
  with arguments
    v[3] - 1-rank numpy array of 3 coordinates representing 3D point
    V[][3] - 2-rank numpy array of vertices  
    T[][3] - 2-rank numpy array of indices of vertices composing triangles
    N[][3] - 2-rank numpy array of normals of triangles
  
  Returns: 
    M[] - 1-rank numpy array of the ratio of the surface that is visible
  
  Ref:
  * http://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.ndarray.html
  * http://docs.scipy.org/doc/numpy/reference/c-api.array.html#creating-arrays
  * http://folk.uio.no/hpl/scripting/doc/python/NumPy/Numeric/numpy-13.html
*/

static PyObject *mesh_visibility(PyObject *self, PyObject *args){


  //
  // Storing/Reading arguments
  //
  
  PyArrayObject *ov, *oV, *oT, *oN;

  // parse arguments
  if (!PyArg_ParseTuple(args, "O!O!O!O!", 
        &PyArray_Type, &ov,
        &PyArray_Type, &oV, 
        &PyArray_Type, &oT,
        &PyArray_Type, &oN)) {
    return NULL;
  }
  
  double *view = (double*)PyArray_DATA(ov);
 
  int 
    Nv = PyArray_DIM(oV, 0),
    Nt = PyArray_DIM(oT, 0);
    
  T3Dpoint<double> *V_begin = (T3Dpoint<double>*)PyArray_DATA(oV);
  std::vector<T3Dpoint<double> > V(V_begin, V_begin + Nv);
  
  Ttriangle *T_begin = (Ttriangle *)PyArray_DATA(oT);
  std::vector<Ttriangle> T(T_begin, T_begin + Nt);
  
  T3Dpoint<double> *N_begin = (T3Dpoint<double>*)PyArray_DATA(oN); 
  std::vector<T3Dpoint<double> > N(N_begin, N_begin + Nt);
  
  std::vector<double> M;

  //
  //  Calculate visibility
  //
  
  triangle_mesh_visibility(view, V, T, N, M);
   
  //
  // Storing result
  //
  
  npy_intp dims[1] = { Nt };
  
  return PyArray_SimpleNewFromVector<double>(1, dims, NPY_DOUBLE, M.data());
}

/*  define functions in module */
static PyMethodDef Methods[] = {
    
    { "roche_critical_potential", 
      roche_critical_potential,   
      METH_VARARGS, 
      "Determine the critical potentials for given values of q, F, and d."},
      
      
    { "roche_pole", 
      roche_pole,   
      METH_VARARGS, 
      "Determine the height of the pole of generalized Roche lobes for given values of q, F, d and Omega0"},
    
    
    // Some modification in declarations due to use of keywords
    // Ref:
    //  * https://docs.python.org/2.0/ext/parseTupleAndKeywords.html
    { "roche_marching_mesh", 
      (PyCFunction)roche_marching_mesh,   
      METH_VARARGS|METH_KEYWORDS, 
      "Determine the triangular meshing of generalized Roche lobes for "
      "given values of q, F, d and value of the generalized Kopal potential "
      "Omega0. The edge of triangles used in the mesh are approximately delta."},
    
    { "mesh_visibility",
      mesh_visibility,
      METH_VARARGS,
      "Determine the ratio of triangle surfaces that are visible in a triangular mesh."},
      
    {NULL,  NULL, 0, NULL} // terminator record
};

/* module initialization */
PyMODINIT_FUNC initlibphoebe (void)
{
  PyObject *backend = Py_InitModule("libphoebe", Methods);
  
  // Added to handle Numpy arrays
  // Ref: 
  // * http://docs.scipy.org/doc/numpy-1.10.1/user/c-info.how-to-extend.html
  import_array();
  
  if (!backend) return;
}