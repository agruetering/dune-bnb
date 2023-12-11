add_library(BB OBJECT bestfirst.cc
	      depthfirst.cc
	      master.cc
	      breadthfirst.cc)
target_compile_options(BB PUBLIC -g)

 