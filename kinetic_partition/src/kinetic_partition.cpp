#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>
#include "boost/tuple/tuple.hpp"
#include "boost/python/object.hpp"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "kinetic_model.h"
#include "lsd_interface.h"
#include "regularization_angles.h"
#include "regularization_angles_ms.h"
#include "regularization_ordinates.h"
#include "regularization_ordinates_ms.h"
#include "propagation.h"
#include "trace.h"
#include "defs.h"

namespace bpn = boost::python::numpy;
namespace bp =  boost::python;
typedef boost::tuple< std::vector< std::vector< float > >, std::vector< std::vector<uint32_t> > > Custom_tuple;

template <class T>
struct VecvecToArray
{//converts a vector< vector<uint32_t> > to a numpy 2d array
    static PyObject * convert(const std::vector< std::vector<T> > & vecvec)
    {
        npy_intp dims[2];
        dims[0] = vecvec.size();
        dims[1] = vecvec[0].size();
        PyObject * obj;
        if (typeid(T) == typeid(uint8_t))
            obj = PyArray_SimpleNew(2, dims, NPY_UINT8);
        else if (typeid(T) == typeid(float))
            obj = PyArray_SimpleNew(2, dims, NPY_FLOAT32);
        else if (typeid(T) == typeid(uint32_t))
            obj = PyArray_SimpleNew(2, dims, NPY_UINT32);
        void * arr_data = PyArray_DATA((PyArrayObject*)obj);
        std::size_t cell_size = sizeof(T);
        for (std::size_t i = 0; i < dims[0]; i++)
        {
            memcpy((char *)arr_data + i * dims[1] * cell_size, &(vecvec[i][0]), dims[1] * cell_size);
        }
        return obj;
    }
};

struct to_py_tuple
{//converts output to a python tuple
    static PyObject* convert(const Custom_tuple& c_tuple){
        bp::list values;
        //add all c_tuple items to "values" list

        PyObject * vertices_vec_pyo = VecvecToArray<float>::convert(c_tuple.get<0>());
        PyObject * edges_vec_pyo = VecvecToArray<uint32_t>::convert(c_tuple.get<1>());

        values.append(bp::handle<>(bp::borrowed(vertices_vec_pyo)));
        values.append(bp::handle<>(bp::borrowed(edges_vec_pyo)));

        return bp::incref( bp::tuple( values ).ptr() );
    }
};

//partition the input image
PyObject* partition_image
(
        const bp::str img_path,
        const float scale = 0.8,
        const float num_intersection = 1,
        const bool enable_regularise = true,
        const bool verbose = false
)
{
    clock_t begin_step= clock();
    if (verbose)
        std::cout<<"Start kinetic partition parameter initialization ... " << std::endl;
    GDALAllRegister();
    Kinetic_Model* model = new Kinetic_Model();
    model->params->path_input_image = std::string(bp::extract<char const *>(img_path));
    model->params->lsd_scale = scale;
    model->params->prop_K = num_intersection;
    model->params->rega_regp_enabled = enable_regularise;

    model->set_time_string();
    model->set_basename();
    model->reinit();

    if (verbose)
        std::cout<<"    Step 1. Runs LSD ... " << std::endl;
    // Step 1. Runs LSD
    LSD::find_segments(model);

    if (model->params->rega_regp_enabled && model->segments.size()>1)
    {
        if (verbose)
            std::cout<<"    Step 2. Regularizes segments' angles ... " << std::endl;
        // Step 2. Regularizes segments' angles
        Regularization_Angles* m_rega = new Regularization_Angles_Mean_Shift();
        m_rega->regularize(model);
        delete m_rega;

        if (verbose)
            std::cout<<"    Step 3. Regularizes segments' ordinates ... " << std::endl;
        // Step 3. Regularizes segments' ordinates
        Regularization_Ordinates* m_regp = new Regularization_Ordinates_Mean_Shift();
        m_regp->regularize(model);
        delete m_regp;
    }
    model->set_angle_hist_index();

    if (verbose)
        std::cout<<"    Step 4. Propagation ... " << std::endl;
    // Step 4. Propagation
    Propagation::propagate_image_domain(model);

    model->elapsed_time_building_graph = double(clock() - begin_step) / CLOCKS_PER_SEC;
    if (verbose)
        std::cout<<"Partition cost " << model->elapsed_time_building_graph << " (s)."<< std::endl;

    if (verbose)
        std::cout<<"Parsing to python datatypes ... " << std::endl;

    // The image coordinates system: origin (0.0, 0.0) at the top-left corner of pixel (0,0), with y-axis pointing down.
    // We transform our coordinates system to the image coordinates.
    int id_vertex = 0;
    int id_edge = 0;
    Point2d shift(0.5, 0.5);
    std::vector< std::vector< float > > vertices(model->graph->v_size, std::vector< float >(2, 0.f));
    Vertex* v = model->graph->vertices_head;
    while (v != NULL)
    {
        v->id_vertex = id_vertex++;
        Point2d pt = v->pt + shift;
        Point2d pc(jclamp(0, pt.x, model->I.cols), jclamp(0, model->I.rows - pt.y, model->I.rows));
        std::vector< float > v_pc = {float(pc.x), float(pc.y)};
        vertices.at( v->id_vertex) = v_pc;
        v = v->v_next;
    }

    // We list the edges
    std::vector< std::vector< uint32_t > > edges(model->graph->e_size, std::vector< uint32_t >(2, -1));
    Edge* e = model->graph->edges_head;
    while (e != NULL)
    {
        e->id_edge = id_edge++;
        std::vector< uint32_t > ed_i = {uint32_t(e->v1->id_vertex), uint32_t(e->v2->id_vertex)};
        edges.at(e->id_edge) = ed_i;
        e = e->e_next;
    }

    delete model;
    return to_py_tuple::convert(Custom_tuple(vertices, edges));
}

BOOST_PYTHON_MODULE(libkinetic_partition)
{
    _import_array();
    Py_Initialize();
    bpn::initialize();
    bp::to_python_converter< Custom_tuple, to_py_tuple>();
    def("partition_image", partition_image);
}

