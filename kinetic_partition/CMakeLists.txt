### CMake config ###
cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0048 NEW)

### Project ###
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

include_directories( BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
endif(MSVC)

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

### External libs ###
add_definitions(-DCGAL_HAS_THREADS)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Python
find_package(Python COMPONENTS Interpreter Development NumPy)

# Boost, numpy + python version
#./bootstrap.sh --with-python=/usr/bin/python3.10
#./b2 --with-python --with-thread --with-system --with-filesystem --with-regex install

# Boost
find_package(Boost 1.65.0 REQUIRED COMPONENTS filesystem numpy)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREAD OFF)
include_directories(${Boost_INCLUDE_DIRS})

# GDAL
find_package(GDAL REQUIRED)
include_directories( BEFORE ${GDAL_INCLUDE_DIRS})

# OpenCV
find_package( OpenCV REQUIRED )

# CGAL
find_package(CGAL COMPONENTS Core REQUIRED)
include( ${CGAL_USE_FILE} )

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

### Complier flags ###
add_definitions( -Wall )

### Source files listing variables ###
set(${PROJECT_NAME}_HEADERS src/defs.h src/comparators.h src/geometry.h
        src/indexed_event.h src/kinetic_model.h src/line_item.h
        src/lsd.h src/lsd_interface.h src/lsd_tuned.h
        src/matrix.h src/means.h src/parameters.h src/partition.h
        src/partition_elements.h src/propagation.h
        src/quadtree.h src/quadtree_point.h src/r_angle.h
        src/refine_partition.h src/regularization_angles.h
        src/regularization_angles_ms.h
        src/regularization_ordinates.h src/regularization_ordinates_ms.h
        src/r_ordinate.h
        src/segment_ray.h src/segment_tree.h src/svg.h src/trace.h)

set( ${PROJECT_NAME}_SOURCES src/geometry.cpp
        src/indexed_event.cpp src/kinetic_model.cpp src/line_item.cpp
        src/lsd.cpp src/lsd_interface.cpp src/lsd_tuned.cpp
        src/means.cpp src/parameters.cpp src/partition.cpp
        src/partition_elements.cpp src/propagation.cpp
        src/quadtree_point.cpp src/r_angle.cpp
        src/refine_partition.cpp src/regularization_angles.cpp
        src/regularization_angles_ms.cpp
        src/regularization_ordinates.cpp src/regularization_ordinates_ms.cpp
        src/r_ordinate.cpp
        src/segment_ray.cpp src/segment_tree.cpp src/kinetic_partition.cpp)

if( CGAL_FOUND AND GDAL_FOUND AND Boost_FOUND AND Eigen3_FOUND AND OpenCV_FOUND )
    ### Library ###
    set(CMAKE_LD_FLAG "${CMAKE_LD_FLAGS} -shared -Wl -fPIC --export-dynamic -o -O3 -Wall -std=c++11")
    add_library(kinetic_partition SHARED ${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_SOURCES} )
    target_link_libraries(kinetic_partition ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
            ${OpenCV_LIBS} ${GDAL_LIBRARIES} Python::Python Python::NumPy ${Boost_LIBRARIES} Boost::filesystem Boost::python Boost::numpy)
    set_target_properties(kinetic_partition PROPERTIES VERSION 1.0.0 SOVERSION 1)
    ### Error messages ###
else()
    message(STATUS "NOTICE: This demo requires ?? and will not be compiled.")
endif()


##########################archived#######################333
#python old
#set(PYTHON_EXECUTABLE /home/geo3d/anaconda3/bin/python3.9)
#set(PYTHON_LIBRARY /home/geo3d/anaconda3/lib/libpython3.9.so)
#set(PYTHON_INCLUDE_DIR /home/geo3d/anaconda3/include/python3.9)
#set(PYTHON_NUMPY_INCLUDE_DIRS /home/geo3d/anaconda3/lib/python3.9/site-packages/numpy/core/include)
#set(PYTHON_NUMPY_VERSION 1.21.5)
#
#set(CMAKE_PREFIX_PATH /home/geo3d/anaconda3/envs/img2roof ${CMAKE_PREFIX_PATH})
#set (conda_env_path "/home/geo3d/anaconda3/")
#list(APPEND CMAKE_PREFIX_PATH ${conda_env_path})

#find_package(PythonLibs)
#find_package(PythonInterp)
#find_package(NumPy REQUIRED)


#set(PYTHONVERSION ${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
#message(STATUS "Found python ${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
#if(Boost_VERSION LESS_EQUAL 106900)
#    message(STATUS "Found Boost ${Boost_VERSION} and BoostNumpy ${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
#else()
#    message(FATAL_ERROR "Boost version ${Boost_VERSION} is not supported. Please install Boost >=1.65.0 and <= 1.69.0.")
#endif()
