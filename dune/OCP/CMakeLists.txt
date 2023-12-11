#install headers
install(FILES probleminterface.hh 
	      heatdriver.hh 
	      adjointdriver.hh
	      linearheatfem.hh
	      adjointheatfem.hh
	      hfunctions.hh
	      switchconstr.hh
	      Dmax.hh
	      ADMMsolver.hh
	      ADMMsystem.hh
 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/OCP)
 

