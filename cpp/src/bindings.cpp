#include "solar_dynamo/model.hpp"
#include <cstring>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl/filesystem.h>

namespace py=pybind11;using namespace solar_dynamo;
namespace{
py::array_t<double> one(const std::vector<double>&v){py::array_t<double>a(v.size());std::memcpy(a.mutable_data(),v.data(),v.size()*sizeof(double));return a;}
py::array_t<double> field(const std::vector<double>&v){py::array_t<double>a({Model::nmax,Model::nmax});std::memcpy(a.mutable_data(),v.data(),v.size()*sizeof(double));return a;}
}
PYBIND11_MODULE(_core,m){
 m.doc()="Step-for-step C++ port of surya_dynamo_v3.f";
 py::class_<Parameters>(m,"Parameters").def(py::init<>())
  .def_readwrite("tmax",&Parameters::tmax).def_readwrite("v0",&Parameters::v0)
  .def_readwrite("et0",&Parameters::et0).def_readwrite("et1",&Parameters::et1)
  .def_readwrite("al0",&Parameters::al0).def_readwrite("omega0",&Parameters::omega0)
  .def_readwrite("dt",&Parameters::dt).def_readwrite("relaxed_initial_state",&Parameters::relaxed_initial_state);
 py::class_<Snapshot>(m,"Snapshot")
  .def_property_readonly("time",[](const Snapshot&s){return s.time;})
  .def_property_readonly("step",[](const Snapshot&s){return s.step;})
  .def_property_readonly("radius",[](const Snapshot&s){return one(s.radius);})
  .def_property_readonly("theta",[](const Snapshot&s){return one(s.theta);})
  .def_property_readonly("poloidal_u",[](const Snapshot&s){return field(s.poloidal_u);})
  .def_property_readonly("poloidal",[](const Snapshot&s){return field(s.poloidal_output);})
  .def_property_readonly("toroidal",[](const Snapshot&s){return field(s.toroidal);})
  .def_property_readonly("omega",[](const Snapshot&s){return field(s.omega);});
 py::class_<Model>(m,"Model").def(py::init<Parameters>(),py::arg("parameters")=Parameters{})
  .def("initialize",&Model::initialize,py::arg("init_file")="init.dat",py::call_guard<py::gil_scoped_release>())
  .def("step",&Model::step,py::arg("count")=1,py::call_guard<py::gil_scoped_release>())
  .def("step_with_output",&Model::step_with_output,py::arg("count"),py::arg("output_directory")=".",py::call_guard<py::gil_scoped_release>())
  .def("run",&Model::run,py::arg("output_directory")=".",py::call_guard<py::gil_scoped_release>())
  .def("write_final_outputs",&Model::write_final_outputs,py::arg("output_directory")=".")
  .def("snapshot",&Model::snapshot).def_property_readonly("time",&Model::time)
  .def_property_readonly("step_number",&Model::step_number).def_property_readonly("initialized",&Model::initialized);
 m.attr("NMAX")=Model::nmax;m.attr("LMAX")=Model::lmax;
}
