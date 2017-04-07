cmake_minimum_required(VERSION 3.8)

#project(my_svd)

include_directories(BEFORE alglib)



#add_executable(main main.cpp)


set (PROJECT my_svd) 
 
project (${PROJECT}) 
 
set (HEADERS 
     alglib_svd.h
     matrix.h
 	 my_bidiag.h
 	 my_bidiag_svd.h
 	 my_svd.h
 	 my_timer.h
 	 my_tridiag_eig.h)
 
set (SOURCES 
 	 alglib_svd.cpp
 	 main.cpp
 	 matrix.cpp
 	 my_bidiag.cpp
 	 my_bidiag_svd.cpp
 	 my_svd.cpp
 	 my_tridiag_eig.cpp)

set (ALGLIB_SRC
	alglib/alglibinternal.cpp
	alglib/alglibmisc.cpp
	alglib/ap.cpp
	alglib/linalg.cpp)

list (APPEND SOURCES ${ALGLIB_SRC})	
 
add_executable (${PROJECT} ${HEADERS} ${SOURCES})

